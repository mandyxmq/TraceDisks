#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <tuple>       // For std::tuple and std::tie
#include <algorithm>   // For std::clamp
#include <omp.h>
#include <chrono>  // Include the chrono library
#include <iomanip>
#include <sstream>
#include <atomic>

// Utility struct for 2D vectors
struct Vec2 {
    double x, y;
    
    Vec2 rotate(double angle) const {
        double cos_theta = std::cos(angle);
        double sin_theta = std::sin(angle);
        return {
            x * cos_theta - y * sin_theta, 
            x * sin_theta + y * cos_theta
        };
    }

    Vec2 operator+(const Vec2& other) const {
        return {x + other.x, y + other.y};
    }

    Vec2 operator-(const Vec2& other) const {
        return {x - other.x, y - other.y};
    }

    Vec2 operator*(double scalar) const {
        return {x * scalar, y * scalar};
    }

    double dot(const Vec2& other) const {
        return x * other.x + y * other.y;
    }

    double norm() const {
        return std::sqrt(x * x + y * y);
    }

    Vec2 normalize() const {
        double n = norm();
        if (n == 0) return {0, 0};        
        return {x / n, y / n};
    }
};

struct Circle {
    Vec2 center;
    double radius;
};

std::mt19937 gen;  // Global random number generator

double random_double(double min, double max) {
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

// std::vector<Circle> create_random_circles(double radius_big, double radius_small, int num_circles) {
//     std::vector<Circle> circles;
//     int attempts = 0;
//     int max_attempts = 1000;
//     double epsilon = 0.1;

//     while (circles.size() < num_circles && attempts < max_attempts) {
//         double angle = random_double(0, 2 * M_PI);
//         double r = random_double(0, radius_big - radius_small);
//         double x = r * cos(angle);
//         double y = r * sin(angle);

//         bool overlap = false;
//         for (const auto& circle : circles) {
//             double distance = std::sqrt((x - circle.center.x) * (x - circle.center.x) +
//                                         (y - circle.center.y) * (y - circle.center.y));
//             if (distance < 2 * radius_small + epsilon) {
//                 overlap = true;
//                 break;
//             }
//         }

//         if (!overlap) {
//             circles.push_back({{x, y}, radius_small});
//         }
//         attempts++;
//     }
//     return circles;
// }

// Function to create random circles using grid-based optimization
std::vector<Circle> create_random_circles(double radius_big, double radius_small, int num_circles) {
    std::vector<Circle> circles;
    int attempts = 0;
    int max_attempts = 1000;
    double epsilon = 0.1;

    // Grid size based on the big circle's radius and the small circles' radius
    int grid_size = static_cast<int>(radius_big / (2 * radius_small));
    std::vector<std::vector<std::vector<Circle>>> grid(grid_size, std::vector<std::vector<Circle>>(grid_size));

    // Helper function to get the grid index for a circle's position
    auto get_grid_index = [&](double x, double y) {
        int ix = std::max(0, std::min(grid_size - 1, static_cast<int>((x + radius_big) / (2 * radius_small))));
        int iy = std::max(0, std::min(grid_size - 1, static_cast<int>((y + radius_big) / (2 * radius_small))));
        return std::make_pair(ix, iy);
    };

    while (circles.size() < num_circles && attempts < max_attempts) {
        // Generate random polar coordinates and convert to Cartesian
        double angle = random_double(0, 2 * M_PI);
        double r = random_double(0, radius_big - radius_small);
        Vec2 center = {r * std::cos(angle), r * std::sin(angle)};

        // Get the grid cell for the new circle
        auto [ix, iy] = get_grid_index(center.x, center.y);

        // Check for overlap with circles in the neighboring cells
        bool overlap = false;
        for (int i = std::max(0, ix - 1); i <= std::min(grid_size - 1, ix + 1); ++i) {
            for (int j = std::max(0, iy - 1); j <= std::min(grid_size - 1, iy + 1); ++j) {
                for (const auto& circle : grid[i][j]) {
                    double distance = (center - circle.center).norm();
                    if (distance < 2 * radius_small + epsilon) {
                        overlap = true;
                        break;
                    }
                }
                if (overlap) break;
            }
            if (overlap) break;
        }

        // If no overlap, add the new circle
        if (!overlap) {
            Circle new_circle = {center, radius_small};
            circles.push_back(new_circle);
            grid[ix][iy].push_back(new_circle);  // Add circle to the grid
        }

        attempts++;
    }

    return circles;
}


Vec2 random_point_on_circle(double radius) {
    double angle = random_double(0, 2 * M_PI);
    return {radius * cos(angle), radius * sin(angle)};
}

double fresnel_reflectance(double cos_i, double n1, double n2) {
    cos_i = std::clamp(cos_i, -1.0, 1.0);
    bool entering = cos_i > 0;
    if (!entering) {
        std::swap(n1, n2);
        cos_i = std::abs(cos_i);
    }

    double sin_t2 = (n1 / n2) * (n1 / n2) * (1 - cos_i * cos_i);
    if (sin_t2 > 1.0 || sin_t2 < 0.0) {
        return 1.0; // Total internal reflection
    }
    //std::cout<<"sin_i*n1: "<<std::sqrt(1 - cos_i * cos_i)*n1<<" sin_t*n2: "<<std::sqrt(sin_t2)*n2<<std::endl;

    double cos_t = std::sqrt(1 - sin_t2);
    double rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t));
    double rp = ((n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t));
    return (rs * rs + rp * rp) / 2;
}

Vec2 reflect(const Vec2& ray_direction, const Vec2& normal) {
    return ray_direction - normal * (2 * ray_direction.dot(normal));
}

Vec2 refract(const Vec2& ray_direction, const Vec2& normal, double n1, double n2) {
    double cos_i = -normal.dot(ray_direction);
    bool entering = cos_i > 0;
    if (!entering) {
        std::swap(n1, n2);
        cos_i = std::abs(cos_i);
    }

    double sin_t2 = (n1 / n2) * (n1 / n2) * (1 - cos_i * cos_i);
    if (sin_t2 > 1.0 || sin_t2 < 0.0) {
        return {0, 0}; // Total internal reflection
    }

    double cos_t = std::sqrt(1 - sin_t2);
    double sin_t = std::sqrt(sin_t2);
    Vec2 perp_normal = ray_direction - normal * ray_direction.dot(normal);
    perp_normal = perp_normal.normalize(); // Normalize to make it a unit vector
    
    Vec2 refract_dir;
    if (entering) {
        refract_dir.x = sin_t * perp_normal.x + cos_t * (-normal.x);
        refract_dir.y = sin_t * perp_normal.y + cos_t * (-normal.y);
    } else {
        refract_dir.x = sin_t * perp_normal.x + cos_t * normal.x;
        refract_dir.y = sin_t * perp_normal.y + cos_t * normal.y;
    } 

    return refract_dir.normalize();
}

bool intersect_ray_circle(const Vec2& ray_origin, const Vec2& ray_direction, const Vec2& circle_center, double circle_radius, Vec2& intersection_point) {
    Vec2 oc = ray_origin - circle_center;
    double a = ray_direction.dot(ray_direction);
    double b = 2.0 * oc.dot(ray_direction);
    double c = oc.dot(oc) - circle_radius * circle_radius;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return false; // No intersection
    }

    double t1 = (-b - std::sqrt(discriminant)) / (2 * a);
    double t2 = (-b + std::sqrt(discriminant)) / (2 * a);
    double t = (t1 > 0) ? t1 : t2;

    if (t < 0) {
        return false; // Intersection is behind the ray origin
    }

    intersection_point = ray_origin + ray_direction * t;
    return true;
}

bool intersect_ray_circle(const Vec2& ray_origin, const Vec2& ray_direction, const Vec2& circle_center, double circle_radius, Vec2& intersection_point, double& t) {
    Vec2 oc = ray_origin - circle_center;
    double a = ray_direction.dot(ray_direction);
    double b = 2.0 * oc.dot(ray_direction);
    double c = oc.dot(oc) - circle_radius * circle_radius;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return false; // No intersection
    }

    double t1 = (-b - std::sqrt(discriminant)) / (2 * a);
    double t2 = (-b + std::sqrt(discriminant)) / (2 * a);
    t = (t1 > 0) ? t1 : t2;

    if (t < 0) {
        return false; // Intersection is behind the ray origin
    }

    intersection_point = ray_origin + ray_direction * t;
    return true;
}

void compute_angular_energy_distribution(const std::vector<double>& outgoing_angles, const std::vector<double>& energies, int num_bins, const std::string& output_file,
        double total_energy) {
    // Initialize bins for angles (-pi to pi) and energies
    std::vector<double> bin_edges(num_bins + 1);
    std::vector<double> bin_centers(num_bins);
    double bin_width = (2 * M_PI) / num_bins;

    // Create bin edges and centers
    for (int i = 0; i <= num_bins; ++i) {
        bin_edges[i] = -M_PI + i * bin_width;
    }
    for (int i = 0; i < num_bins; ++i) {
        bin_centers[i] = (bin_edges[i] + bin_edges[i + 1]) / 2;
    }

    // Accumulate energy into bins
    std::vector<double> accumulated_energy(num_bins, 0.0);

    for (size_t i = 0; i < outgoing_angles.size(); ++i) {
        // Adjust angle to the range (-pi, pi]
        double angle = std::fmod(outgoing_angles[i] + 2 * M_PI, 2 * M_PI) - M_PI;

        // Find the appropriate bin for the current angle
        int bin_index = static_cast<int>((angle + M_PI) / bin_width);
        if (bin_index >= 0 && bin_index < num_bins) {
            accumulated_energy[bin_index] += energies[i] / total_energy;
        }
    }

    // Save the accumulated energy to a binary file
    std::ofstream outfile(output_file, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << output_file << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<char*>(accumulated_energy.data()), num_bins * sizeof(double));
    outfile.close();
}

void compute_distributions(const std::vector<double>& outgoing_angles, const std::vector<double>& energies, int num_bins, std::vector<double>& accumulated, double total_energy) {
    // Initialize bins for angles (-pi to pi) and energies
    std::vector<double> bin_edges(num_bins + 1);
    std::vector<double> bin_centers(num_bins);
    double bin_width = (2 * M_PI) / num_bins;

    // Create bin edges and centers
    for (int i = 0; i <= num_bins; ++i) {
        bin_edges[i] = -M_PI + i * bin_width;
    }
    for (int i = 0; i < num_bins; ++i) {
        bin_centers[i] = (bin_edges[i] + bin_edges[i + 1]) / 2;
    }

    for (size_t i = 0; i < outgoing_angles.size(); ++i) {
        // Adjust angle to the range (-pi, pi]
        double angle = std::fmod(outgoing_angles[i] + 2 * M_PI, 2 * M_PI) - M_PI;

        // Find the appropriate bin for the current angle
        int bin_index = static_cast<int>((angle + M_PI) / bin_width);
        if (bin_index >= 0 && bin_index < num_bins) {
            accumulated[bin_index] += energies[i];
        }
    }
}

void compute_distributions(const std::vector<double>& outgoing_angles, const std::vector<double>& energies, int num_bins, std::vector<double>& accumulated) {
    // Initialize bins for angles (-pi to pi) and energies
    std::vector<double> bin_edges(num_bins + 1);
    std::vector<double> bin_centers(num_bins);
    double bin_width = (2 * M_PI) / num_bins;

    std::vector<double> count(num_bins, 0.0);

    // Create bin edges and centers
    for (int i = 0; i <= num_bins; ++i) {
        bin_edges[i] = -M_PI + i * bin_width;
    }
    for (int i = 0; i < num_bins; ++i) {
        bin_centers[i] = (bin_edges[i] + bin_edges[i + 1]) / 2;
    }

    for (size_t i = 0; i < outgoing_angles.size(); ++i) {
        // Adjust angle to the range (-pi, pi]
        double angle = std::fmod(outgoing_angles[i] + 2 * M_PI, 2 * M_PI) - M_PI;

        // Find the appropriate bin for the current angle
        int bin_index = static_cast<int>((angle + M_PI) / bin_width);
        if (bin_index >= 0 && bin_index < num_bins) {
            accumulated[bin_index] += energies[i];
            count[bin_index] += 1.0;
        }
    }

    // normalize by count
    for (int i = 0; i < num_bins; ++i) {
        if (count[i] > 0) {
            accumulated[i] /= count[i];
        }
    }
}


void trace_ray(const Vec2& start_point, Vec2 ray_direction, 
                const std::vector<Circle>& circles, double radius_big, double n1, double n2, double sigma_a, double incident_angle, 
                std::vector<std::tuple<double, double, Vec2>>& results, std::vector<Vec2>& interaction_points) {
    double energy = 1.0;
    Vec2 ray_origin = start_point;
    bool interacted = false;
    double lowest_energy = 1e-4;
    double epsilon = 1e-6;

    for (int i = 0; i < 100; ++i) {
        Vec2 closest_intersection;
        Circle closest_circle = {{0, 0}, 0};
        Vec2 closest_normal;
        bool intersection_found = false;

        // prevent self-intersection
        ray_origin = ray_origin + ray_direction * epsilon;

        for (const auto& circle : circles) {
            Vec2 intersection;
            if (intersect_ray_circle(ray_origin, ray_direction, circle.center, circle.radius, intersection)) {
                double distance = (intersection - ray_origin).norm();
                if (!intersection_found || distance < (closest_intersection - ray_origin).norm()) {
                    closest_intersection = intersection;
                    closest_circle = circle;
                    closest_normal = (intersection - circle.center).normalize();
                    intersection_found = true;
                }
            }
        }

        if (!intersection_found) {
            break;
        }

        interacted = true;
        double entry_to_intersection_length = (closest_intersection - ray_origin).norm();
        double transmission = std::exp(-sigma_a * entry_to_intersection_length / (2 * closest_circle.radius));
        energy *= transmission;

        if (energy < lowest_energy) {
            return; // Path terminated with negligible energy
        }

        double cos_i = -ray_direction.dot(closest_normal);
        double fresnel = fresnel_reflectance(cos_i, n1, n2);
        double random_number = random_double(0, 1);
        //std::cout<<"fresnel: "<<fresnel<<", random_number: "<<random_number<<std::endl;
        //std::cout<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;

        if (random_number < fresnel) {
            ray_direction = reflect(ray_direction, closest_normal);
        } else {
            Vec2 refracted_direction = refract(ray_direction, closest_normal, n1, n2);
            if (refracted_direction.norm() == 0) {
                ray_direction = reflect(ray_direction, closest_normal); // Total internal reflection
            } else {
                ray_direction = refracted_direction;
            }
        }

        ray_origin = closest_intersection;

        // if (interacted)
        //     // Store the interaction point
        //     interaction_points.push_back(closest_intersection);
    }

    if (interacted) {
        double outgoing_angle = std::atan2(ray_direction.y, ray_direction.x) - incident_angle;
        outgoing_angle = fmod(outgoing_angle + 2 * M_PI, 2 * M_PI);
        results.push_back({outgoing_angle, energy, ray_origin});
        
    }
}

bool trace_ray(const Vec2& start_point, Vec2 ray_direction, 
                const std::vector<Circle>& circles, double radius_big, double n1, double n2, double sigma_a, double incident_angle, 
                std::vector<double>& outgoing_angles, std::vector<double>& energies,
                std::vector<double>& path_lengths, std::vector<double>& Rs, std::vector<double>& Ts) {
    double energy = 1.0;
    Vec2 ray_origin = start_point;
    bool interacted = false;
    double lowest_energy = 1e-4;
    double epsilon = 1e-6;
    double path_length = 0.0;
    double R = 0.0;
    double T = 0.0;

    for (int i = 0; i < 100; ++i) {
        Vec2 closest_intersection;
        Circle closest_circle = {{0, 0}, 0};
        Vec2 closest_normal;
        bool intersection_found = false;

        // prevent self-intersection
        ray_origin = ray_origin + ray_direction * epsilon;

        for (const auto& circle : circles) {
            Vec2 intersection;
            if (intersect_ray_circle(ray_origin, ray_direction, circle.center, circle.radius, intersection)) {
                double distance = (intersection - ray_origin).norm();
                if (!intersection_found || distance < (closest_intersection - ray_origin).norm()) {
                    closest_intersection = intersection;
                    closest_circle = circle;
                    closest_normal = (intersection - circle.center).normalize();
                    intersection_found = true;
                }
            }
        }

        if (!intersection_found) {
            break;
        }

        interacted = true;
        double entry_to_intersection_length = (closest_intersection - ray_origin).norm();
        double transmission = std::exp(-sigma_a * entry_to_intersection_length / (2 * radius_big));
        energy *= transmission;

        if (energy < lowest_energy) {
            return interacted; // Path terminated with negligible energy
        }

        double cos_i = -ray_direction.dot(closest_normal);
        double fresnel = fresnel_reflectance(cos_i, n1, n2);
        double random_number = random_double(0, 1);
        //std::cout<<"fresnel: "<<fresnel<<", random_number: "<<random_number<<std::endl;
        //std::cout<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;

        if (random_number < fresnel) {
            ray_direction = reflect(ray_direction, closest_normal);
            R += 1.0;
        } else {
            Vec2 refracted_direction = refract(ray_direction, closest_normal, n1, n2);
            if (refracted_direction.norm() == 0) {
                ray_direction = reflect(ray_direction, closest_normal); // Total internal reflection
            } else {
                ray_direction = refracted_direction;
            }
            T += 1.0;
        }

        // compute path length
        path_length += entry_to_intersection_length;

        ray_origin = closest_intersection;

        // if (interacted)
        //     // Store the interaction point
        //     interaction_points.push_back(closest_intersection);
    }

    if (interacted) {
        double outgoing_angle = std::atan2(ray_direction.y, ray_direction.x) - incident_angle;
        outgoing_angle = fmod(outgoing_angle + 2 * M_PI, 2 * M_PI);
        //results.push_back({outgoing_angle, energy, ray_origin});
        outgoing_angles.push_back(outgoing_angle);
        energies.push_back(energy);
        path_lengths.push_back(path_length);
        Rs.push_back(R);
        Ts.push_back(T);        
    }

    return interacted;
}

// Function to rotate all circles around the origin by a given angle
void rotate_circles(std::vector<Circle>& circles, double angle) {
    for (auto& circle : circles) {
        circle.center = circle.center.rotate(angle);
    }
}

// void simulate_small_circles(int num_paths, double radius_big, double radius_small, 
//             int num_circles, double n1, double n2, double sigma_a, unsigned int seed, int num_bins) {
//     // Set the random seed
//     gen.seed(seed);

//     std::vector<Circle> circles = create_random_circles(radius_big, radius_small, num_circles);
    
//     // // Store origins of small circles
//     // std::ofstream circle_file("circle_origins.csv");
//     // circle_file << "Center X,Center Y\n";
//     // for (const auto& circle : circles) {
//     //     circle_file << circle.center.x << "," << circle.center.y << "\n";
//     // }
//     // circle_file.close();

//     //std::vector<Vec2> interaction_points;
//     //Vec2 start_point = random_point_on_circle(radius_big);
//     // Vec2 start_point = {-10, 0};
//     //double incident_angle = random_double(0, 2 * M_PI) + M_PI;
//     double incident_angle = 0.0;
//     //std::vector<std::tuple<double, double, Vec2>> results;
//     std::vector<double> energies;
//     std::vector<double> outgoing_angles;

//     double epsilon = M_PI / 100.0;
//     #pragma omp parallel for
//     for (int i = 0; i < num_paths; ++i) {

//         double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_paths));  
//         //double y_start = -8.0;
//         // Uniformly distribute starting points
//         Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Start just left of the circle
//         Vec2 ray_direction = {std::cos(incident_angle), std::sin(incident_angle)};
        
//         double energy = 1.0;  // Initial energy
//         bool has_intersected = false;
//         bool exit = false;
//         int j = 0;

//         double theta = (M_PI - 2 * epsilon) / num_paths * i + M_PI / 2 + epsilon;
//         Vec2 start_point = {std::cos(theta)*radius_big, std::sin(theta)*radius_big};
//         trace_ray(ray_origin, ray_direction, circles, radius_big, n1, n2, 
//                     sigma_a, incident_angle, outgoing_angles, energies);
//         // #pragma omp critical
//         // {
//         //     interaction_points.insert(interaction_points.end(), local_interaction_points.begin(), local_interaction_points.end());
//         // }
//     }

//     // // Write interaction points to a file
//     // std::ofstream interaction_file("interaction_points.csv");
//     // interaction_file << "X,Y\n";
//     // for (const auto& point : interaction_points) {
//     //     interaction_file << point.x << "," << point.y << "\n";
//     // }
//     // interaction_file.close();

//     // std::ofstream outfile("ray_simulation_results.csv");
//     // outfile << "Outgoing Angle (radians),Energy,Exit Point X,Exit Point Y\n";
//     // for (const auto& result : results) {
//     //     double angle, energy;
//     //     Vec2 exit_point;
//     //     std::tie(angle, energy, exit_point) = result;
//     //     outfile << angle << "," << energy << "," << exit_point.x << "," << exit_point.y << "\n";
//     // }
//     // outfile.close();

//      // Save the energy distribution to a binary file    
//     compute_angular_energy_distribution(outgoing_angles, energies, num_bins, "multiple.bin");
// }


void simulate_small_circles(int num_paths, double radius_big, double radius_small, 
                            int num_circles, double n1, double n2, double sigma_a, 
                            unsigned int seed, int num_bins, int num_exp,
                            std::string output_file) {
    // Set the random seed (use thread-private random generators later)
    gen.seed(seed);

    //double total_energy = num_paths * num_exp;
    double total_energy = 0.0;

    // std::vector<Circle> circles = create_random_circles(radius_big, radius_small, num_circles);

    double incident_angle = 0.0;
    double epsilon = M_PI / 100.0;

    // rotate circles and trace rays
    double angle = 2 * M_PI / num_exp;    
    std::vector<double> accumulated_energy(num_bins, 0.0);
    std::vector<double> accumulated_path_length(num_bins, 0.0);
    std::vector<double> accumulated_R(num_bins, 0.0);
    std::vector<double> accumulated_T(num_bins, 0.0);
    
    for (int index = 0; index < num_exp; ++index){
        std::cout<<"exp: "<<index<<std::endl;
        // rotate circles
        // rotate_circles(circles, angle);

        std::vector<Circle> circles = create_random_circles(radius_big, radius_small, num_circles);

        // Shared vectors to accumulate results
        std::vector<double> outgoing_angles;
        std::vector<double> energies;
        std::vector<double> path_lengths;
        std::vector<double> Rs;
        std::vector<double> Ts;

        // Use parallel region with private local vectors
        #pragma omp parallel
        {
            // Thread-private vectors
            std::vector<double> local_outgoing_angles;
            std::vector<double> local_energies;
            std::vector<double> local_path_length;
            std::vector<double> local_R;
            std::vector<double> local_T;

            #pragma omp for
            for (int i = 0; i < num_paths; ++i) {
                double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_paths));  
                Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Start just left of the circle
                Vec2 ray_direction = {std::cos(incident_angle), std::sin(incident_angle)};
                
                double theta = (M_PI - 2 * epsilon) / num_paths * i + M_PI / 2 + epsilon;
                Vec2 start_point = {std::cos(theta) * radius_big, std::sin(theta) * radius_big};
                
                // // Simulate the ray and store results in local vectors
                // trace_ray(ray_origin, ray_direction, circles, radius_big, n1, n2, 
                //         sigma_a, incident_angle, local_outgoing_angles, local_energies,
                //         local_path_length, local_R, local_T);
                
                bool interacted = trace_ray(ray_origin, ray_direction, circles, radius_big, n1, n2, sigma_a, 
                              incident_angle, local_outgoing_angles, local_energies, 
                              local_path_length, local_R, local_T);
               
                if (interacted) {
                    #pragma omp critical
                    {
                        total_energy += 1.0;
                    }
                }
            }

            // Merge local results into shared vectors using critical sections
            #pragma omp critical
            {
                outgoing_angles.insert(outgoing_angles.end(), local_outgoing_angles.begin(), local_outgoing_angles.end());
                energies.insert(energies.end(), local_energies.begin(), local_energies.end());
                path_lengths.insert(path_lengths.end(), local_path_length.begin(), local_path_length.end());
                Rs.insert(Rs.end(), local_R.begin(), local_R.end());
                Ts.insert(Ts.end(), local_T.begin(), local_T.end());
            }
        }

        // // Save the energy distribution to a binary file
        // std::stringstream ss;
        // // Set the width to 3, fill with '0', and convert the number to a string
        // ss << std::setw(3) << std::setfill('0') << index;
        // std::string index_str = ss.str();
        // // std::cout<<"index_str: "<<index_str<<std::endl;
        // //std::string index_str = std::to_string(index);
        // compute_angular_energy_distribution(outgoing_angles, energies, num_bins, "multiple_"+index_str+".bin");

        // compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy, total_energy);
        // compute_distributions(outgoing_angles, path_lengths, num_bins, accumulated_path_length);
        // compute_distributions(outgoing_angles, Rs, num_bins, accumulated_R);
        // compute_distributions(outgoing_angles, Ts, num_bins, accumulated_T);
        compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy, total_energy);
        compute_distributions(outgoing_angles, path_lengths, num_bins, accumulated_path_length);
        compute_distributions(outgoing_angles, Rs, num_bins, accumulated_R);
        compute_distributions(outgoing_angles, Ts, num_bins, accumulated_T);
    }

    //std::cout<<"total_energy: "<<total_energy<<" num_paths*num_exp "<< num_paths*num_exp << std::endl;
    // normalize accumulated_energy
    for (int i = 0; i < num_bins; ++i) {
        accumulated_energy[i] /= total_energy;
    }

    // Save the accumulated energy to a binary file
    std::ofstream outfile(output_file, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << output_file << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<char*>(accumulated_energy.data()), num_bins * sizeof(double));
    outfile.close();

    // // save path length to a binary file
    // // add path_length to the output_file
    // std::ofstream outfile_path_length("path_length_"+output_file, std::ios::binary);
    // if (!outfile_path_length) {
    //     std::cerr << "Error opening file for writing: " << "path_length_"+output_file << std::endl;
    //     return;
    // }
    // outfile_path_length.write(reinterpret_cast<char*>(accumulated_path_length.data()), num_bins * sizeof(double));
    // outfile_path_length.close();

    // // save R distribution to a binary file
    // // add R to the output_file
    // std::ofstream outfile_R("R_"+output_file, std::ios::binary);
    // if (!outfile_R) {
    //     std::cerr << "Error opening file for writing: " << "R_"+output_file << std::endl;
    //     return;
    // }
    // outfile_R.write(reinterpret_cast<char*>(accumulated_R.data()), num_bins * sizeof(double));
    // outfile_R.close();

    // // save T distribution to a binary file
    // // add T to the output_file
    // std::ofstream outfile_T("T_"+output_file, std::ios::binary);
    // if (!outfile_T) {
    //     std::cerr << "Error opening file for writing: " << "T_"+output_file << std::endl;
    //     return;
    // }
    // outfile_T.write(reinterpret_cast<char*>(accumulated_T.data()), num_bins * sizeof(double));
    // outfile_T.close();
}

void simulate_small_circles_h(int num_paths, double radius_big, double radius_small, 
                            int num_circles, double n1, double n2, double sigma_a, 
                            unsigned int seed, int num_bins, int num_exp,
                            std::string output_file, double h) {
    // Set the random seed (use thread-private random generators later)
    gen.seed(seed);

    //double total_energy = num_paths * num_exp;
    double total_energy = 0.0;

    // std::vector<Circle> circles = create_random_circles(radius_big, radius_small, num_circles);

    double incident_angle = 0.0;
    double epsilon = M_PI / 100.0;

    // rotate circles and trace rays
    double angle = 2 * M_PI / num_exp;    
    std::vector<double> accumulated_energy(num_bins, 0.0);
    std::vector<double> accumulated_path_length(num_bins, 0.0);
    std::vector<double> accumulated_R(num_bins, 0.0);
    std::vector<double> accumulated_T(num_bins, 0.0);
    
    for (int index = 0; index < num_exp; ++index){
        std::cout<<"exp: "<<index<<std::endl;
        // rotate circles
        // rotate_circles(circles, angle);

        std::vector<Circle> circles = create_random_circles(radius_big, radius_small, num_circles);

        // Shared vectors to accumulate results
        std::vector<double> outgoing_angles;
        std::vector<double> energies;
        std::vector<double> path_lengths;
        std::vector<double> Rs;
        std::vector<double> Ts;

        // Use parallel region with private local vectors
        #pragma omp parallel
        {
            // Thread-private vectors
            std::vector<double> local_outgoing_angles;
            std::vector<double> local_energies;
            std::vector<double> local_path_length;
            std::vector<double> local_R;
            std::vector<double> local_T;

            double y_start = -radius_big + (2.0 * radius_big * (h + 1) / 2.0);  
            Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Start just left of the circle
            Vec2 ray_direction = {std::cos(incident_angle), std::sin(incident_angle)};
            
            // // Simulate the ray and store results in local vectors
            // trace_ray(ray_origin, ray_direction, circles, radius_big, n1, n2, 
            //         sigma_a, incident_angle, local_outgoing_angles, local_energies,
            //         local_path_length, local_R, local_T);
            
            bool interacted = trace_ray(ray_origin, ray_direction, circles, radius_big, n1, n2, sigma_a, 
                            incident_angle, local_outgoing_angles, local_energies, 
                            local_path_length, local_R, local_T);
            
            if (interacted) {
                #pragma omp critical
                {
                    total_energy += 1.0;
                }
            }

            // Merge local results into shared vectors using critical sections
            #pragma omp critical
            {
                outgoing_angles.insert(outgoing_angles.end(), local_outgoing_angles.begin(), local_outgoing_angles.end());
                energies.insert(energies.end(), local_energies.begin(), local_energies.end());
                path_lengths.insert(path_lengths.end(), local_path_length.begin(), local_path_length.end());
                Rs.insert(Rs.end(), local_R.begin(), local_R.end());
                Ts.insert(Ts.end(), local_T.begin(), local_T.end());
            }
        }

        // // Save the energy distribution to a binary file
        // std::stringstream ss;
        // // Set the width to 3, fill with '0', and convert the number to a string
        // ss << std::setw(3) << std::setfill('0') << index;
        // std::string index_str = ss.str();
        // // std::cout<<"index_str: "<<index_str<<std::endl;
        // //std::string index_str = std::to_string(index);
        // compute_angular_energy_distribution(outgoing_angles, energies, num_bins, "multiple_"+index_str+".bin");

        compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy, total_energy);
        compute_distributions(outgoing_angles, path_lengths, num_bins, accumulated_path_length);
        compute_distributions(outgoing_angles, Rs, num_bins, accumulated_R);
        compute_distributions(outgoing_angles, Ts, num_bins, accumulated_T);
    }

    //std::cout<<"total_energy: "<<total_energy<<" num_paths*num_exp "<< num_paths*num_exp << std::endl;
    // normalize accumulated_energy
    for (int i = 0; i < num_bins; ++i) {
        accumulated_energy[i] /= total_energy;
    }

    // Save the accumulated energy to a binary file
    std::ofstream outfile(output_file, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << output_file << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<char*>(accumulated_energy.data()), num_bins * sizeof(double));
    outfile.close();

    // // save path length to a binary file
    // // add path_length to the output_file
    // std::ofstream outfile_path_length("path_length_"+output_file, std::ios::binary);
    // if (!outfile_path_length) {
    //     std::cerr << "Error opening file for writing: " << "path_length_"+output_file << std::endl;
    //     return;
    // }
    // outfile_path_length.write(reinterpret_cast<char*>(accumulated_path_length.data()), num_bins * sizeof(double));
    // outfile_path_length.close();

    // // save R distribution to a binary file
    // // add R to the output_file
    // std::ofstream outfile_R("R_"+output_file, std::ios::binary);
    // if (!outfile_R) {
    //     std::cerr << "Error opening file for writing: " << "R_"+output_file << std::endl;
    //     return;
    // }
    // outfile_R.write(reinterpret_cast<char*>(accumulated_R.data()), num_bins * sizeof(double));
    // outfile_R.close();

    // // save T distribution to a binary file
    // // add T to the output_file
    // std::ofstream outfile_T("T_"+output_file, std::ios::binary);
    // if (!outfile_T) {
    //     std::cerr << "Error opening file for writing: " << "T_"+output_file << std::endl;
    //     return;
    // }
    // outfile_T.write(reinterpret_cast<char*>(accumulated_T.data()), num_bins * sizeof(double));
    // outfile_T.close();
}

void simulate_large_circle(int num_rays, double radius_big, double n1, double n2, double sigma_a, unsigned int seed, int max_bounces, int num_bins, std::string output_file) {
    std::mt19937 gen(seed);
    Vec2 circle_center = {0.0, 0.0};  // Center of the big circle
    std::vector<std::tuple<double, double, double, Vec2, double>> results;  // Store outgoing angle and energy
    //std::vector<std::tuple<double, double, double, Vec2, Vec2, double>> results;  // Store outgoing angle and energy
    std::vector<double> energies;
    std::vector<double> outgoing_angles;
    double min_energy = 1e-4;
    double epsilon = 1e-6;
    std::vector<double> accumulated_energy(num_bins, 0.0);

    // Emit parallel rays from the left side of the circle
    //#pragma omp parallel for
    for (int i = 0; i < num_rays; ++i) {
        double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_rays)) + 1.0 * radius_big / double(num_rays);  
        //double y_start = -8.0;
        // Uniformly distribute starting points
        Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Start just left of the circle
        Vec2 ray_direction = {1.0, 0.0};  // Horizontal ray moving rightward
        
        Vec2 current_origin = ray_origin;
        Vec2 previous_origin = ray_origin;
        Vec2 current_direction = ray_direction;
        double energy = 1.0;  // Initial energy
        bool has_intersected = false;
        bool exit = false;
        int j = 0;

        for (j = 0; j < max_bounces; ++j) {  
            Vec2 intersection_point;
            if (!intersect_ray_circle(current_origin, current_direction, circle_center, radius_big, intersection_point)) {
                exit = true;
                break;  // No intersection with the circle; ray exits
            }
            
            // std::cout<<"intersection_point: "<<intersection_point.x<<", "<<intersection_point.y<<std::endl;

            has_intersected = true;
            // Calculate normal at the intersection point
            Vec2 normal = (intersection_point - circle_center).normalize();

            // Calculate path length inside the circle for absorption
            if (j>0){
                double path_length = (intersection_point - current_origin).norm();
                double transmission = std::exp(-sigma_a * path_length / (2 * radius_big));
                energy *= transmission;
            }            

            // Check if energy is too low to continue
            if (energy < min_energy) {
                break;  // Terminate path
            }

            // Determine reflection or refraction at the circle's surface
            double cos_i = -ray_direction.dot(normal);
            double fresnel = fresnel_reflectance(cos_i, n1, n2);
            //double fresnel = 1.0;
            double random_number = random_double(0, 1);

            //std::cout<<"fresnel: "<<fresnel<<", random_number: "<<random_number<<std::endl;

            if (random_number < fresnel) {
                // Reflect the ray
                // std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                // std::cout<<"normal: "<<normal.x<<", "<<normal.y<<std::endl;
                current_direction = current_direction - normal * (2 * current_direction.dot(normal));
                //std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                // // compute reflected angle and see if it equals to the incident angle
                // double reflected_angle = current_direction.dot(normal);
                // std::cout<<"incident angle "<<cos_i<<", reflected angle "<<reflected_angle<<std::endl;
                //abort();

                // fresnel = fresnel_reflectance(cos_i, n1, n2);
                // energy *= fresnel;
            } else {
                // Refract the ray
                Vec2 refracted_direction = refract(current_direction, normal, n1, n2);
                // Vec2 tmp_dir;
                // tmp_dir.x = -0.736996, tmp_dir.y = -0.675898;
                // Vec2 tmp_normal;
                // tmp_normal.x = 0.736996, tmp_normal.y = 0.675898;
                // Vec2 refracted_direction = refract(tmp_dir, tmp_normal, 1.0, 1.53);
                // if (std::isnan(refracted_direction.x) || std::isnan(refracted_direction.y)) {
                //     std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                //     std::cout<<"normal: "<<normal.x<<", "<<normal.y<<std::endl;
                //     std::cout<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;
                //     abort();
                // }
                if (refracted_direction.norm() == 0) {
                    // Total internal reflection
                    current_direction = current_direction - normal * (2 * current_direction.dot(normal));
                } else {
                    current_direction = refracted_direction;
                }
            }

            // Update the ray's origin to the intersection point
            previous_origin = current_origin;
            current_origin = intersection_point + current_direction * epsilon;

        }

        
        // Store the outgoing angle and energy if the ray interacted with the circle
        // std::cout<<"has_intersected: "<<has_intersected<<", exit: "<<exit<<std::endl;
        if (has_intersected && exit) {
            double outgoing_angle = std::atan2(current_direction.y, current_direction.x);
            #pragma omp critical
            {
                results.emplace_back(outgoing_angle, energy, y_start, current_origin, double(j));
                energies.push_back(energy);
                outgoing_angles.push_back(outgoing_angle);
            }
        }

        // // Store the outgoing angle and energy if the ray interacted with the circle
        // // std::cout<<"has_intersected: "<<has_intersected<<", exit: "<<exit<<std::endl;
        // if (has_intersected && exit && j==3) {
        //     double outgoing_angle = std::atan2(current_direction.y, current_direction.x);
        //     #pragma omp critical
        //     {
        //         results.emplace_back(outgoing_angle, energy, y_start, previous_origin, current_origin, double(j));
        //         energies.push_back(energy);
        //         outgoing_angles.push_back(outgoing_angle);
        //     }
        // }
    }

    // // Write the results to a CSV file
    // std::ofstream outfile("large_circle_simulation_results.csv");
    // outfile << "Outgoing Angle (radians),Energy,Start Y,Exit Point X,Exit Point Y,Bounce\n";
    // for (const auto& result : results) {
    //     double angle, energy, y_start, bounce;
    //     Vec2 exit_point;
    //     std::tie(angle, energy, y_start, exit_point, bounce) = result;
    //     outfile << angle << "," << energy << "," << y_start << "," << exit_point.x << "," << exit_point.y << "," << bounce << "\n";
    // }
    // outfile.close();

    // // Write the results to a CSV file
    // std::ofstream outfile("large_circle_simulation_results.csv");
    // outfile << "Outgoing Angle (radians),Energy,Start Y,Prev Point X,Prev Point Y,Exit Point X,Exit Point Y,Bounce\n";
    // for (const auto& result : results) {
    //     double angle, energy, y_start, bounce;
    //     Vec2 previous_point;
    //     Vec2 exit_point;
    //     std::tie(angle, energy, y_start, previous_point, exit_point, bounce) = result;
    //     outfile << angle << "," << energy << "," << y_start << "," << previous_point.x <<
    //     "," << previous_point.y << "," << exit_point.x << "," << exit_point.y << "," << bounce << "\n";
    // }
    // outfile.close();

    // Save the energy distribution to a binary file    
    //compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy);
    compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy, num_rays);
    // Save the accumulated energy to a binary file
    std::ofstream outfile(output_file, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << output_file << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<char*>(accumulated_energy.data()), num_bins * sizeof(double));
    outfile.close();
}


void simulate_double_circles(int num_rays, double radius_big, double radius_small, 
                            double n_air, double n_big, double n_small,
                            double sigma_a_big, double sigma_a_small,
                            unsigned int seed, int max_bounces, int num_bins, std::string output_file) {
    std::mt19937 gen(seed);
    Vec2 circle_center = {0.0, 0.0};  // Center of the big circle
    std::vector<std::tuple<double, double, double, Vec2, double>> results;  // Store outgoing angle and energy
    std::vector<double> energies;
    std::vector<double> outgoing_angles;
    double min_energy = 1e-4;
    double epsilon = 1e-6;
    std::vector<double> accumulated_energy(num_bins, 0.0);

    // Emit parallel rays from the left side of the circle
    // #pragma omp parallel for
    for (int i = 0; i < num_rays; ++i) {
        //std::cout<<"ray: "<<i<<std::endl;
        double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_rays)) + 1.0 * radius_big / double(num_rays);  
        //double y_start = 4;
        //std::cout<<"y_start: "<<y_start<<std::endl;
        //double y_start = 0.0;
        // Uniformly distribute starting points
        Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Start just left of the circle
        Vec2 ray_direction = {1.0, 0.0};  // Horizontal ray moving rightward
        
        Vec2 current_origin = ray_origin;
        Vec2 previous_origin = ray_origin;
        Vec2 current_direction = ray_direction;
        double energy = 1.0;  // Initial energy
        bool has_intersected = false;
        bool exit = false;
        int j = 0;

        int medium = 0; // 0: air, 1: big circle, 2: small circle
        double n1 = 0.0, n2 = 0.0;
        bool reflect = true;

        for (j = 0; j < 100; ++j) {  
            Vec2 intersection_point, intersection_point_small, intersection_point_big;
            bool intersected_big = false;
            bool intersected_small = false;
            double t = -1.0, t_large = -1.0, t_small = -1.0; // ray distance
            if (intersect_ray_circle(current_origin, current_direction, circle_center, radius_big, intersection_point_big, t_large))
                intersected_big = true;
            if (intersect_ray_circle(current_origin, current_direction, circle_center, radius_small, intersection_point_small, t_small))
                intersected_small = true;

            // std::cout<<"intersected_big: "<<intersected_big<<", intersected_small: "<<intersected_small<<std::endl;
            // std::cout<<"t_large: "<<t_large<<", t_small: "<<t_small<<std::endl;

            if (!(intersected_big || intersected_small)) {
                exit = true;
                break;  // No intersection with the circle; ray exits
            }               

            //std::cout<<"t_large: "<<t_large<<", t_small: "<<t_small<<std::endl;

            has_intersected = true;
            if (t_large < 0){
                t = t_small;
                intersection_point = intersection_point_small;
                intersected_big = false;
                intersected_small = true;
            }                
            else if (t_small < 0){
                t = t_large;
                intersection_point = intersection_point_big;
                intersected_big = true;
                intersected_small = false;                                
            }                
            else{
                if (t_large < t_small){
                    t = t_large;
                    intersection_point = intersection_point_big;       
                    intersected_big = true;
                    intersected_small = false;
                }                    
                else{
                    t = t_small;
                    intersection_point = intersection_point_small;
                    intersected_big = false;
                    intersected_small = true;
                }
            }


            //std::cout<<"intersection_point "<<intersection_point.x<<", "<<intersection_point.y<<std::endl;
            
            // Calculate normal at the intersection point
            Vec2 normal = (intersection_point - circle_center).normalize();

            // Calculate path length inside the circle for absorption
            if (j>0){
                double path_length = (intersection_point - current_origin).norm();
                // if (in_big_circle) use sigma_a_big, else use sigma_a_small
                double sigma_a = 0.0;
                if (medium == 1)
                    sigma_a = sigma_a_big;
                else if (medium == 2)
                    sigma_a = sigma_a_small;
                else
                    sigma_a = 0.0;
                double transmission = std::exp(-sigma_a * path_length / (2 * radius_big));
                energy *= transmission;
            }            

            // Check if energy is too low to continue
            if (energy < min_energy) {
                break;  // Terminate path
            }

            // Determine reflection or refraction at the circle's surface
            double cos_i = -ray_direction.dot(normal);
            // if in big circle and intersect with big circle, n1 = n_air, n2 = n_big
            // if in big circle and intersect with small circle, n1 = n_big, n2 = n_small
            // if in small circle, n1 = n_big, n2 = n_small

            // std::cout<<"in_big_circle: "<<in_big_circle<<std::endl;
            // std::cout<<"intersected_big: "<<intersected_big<<", intersected_small: "<<intersected_small<<std::endl;

            if (medium == 0){
                n1 = n_air;
                n2 = n_big;
            } else if (medium == 1){
                if (intersected_small){
                    n1 = n_big;
                    n2 = n_small;
                } else{
                    n1 = n_air;
                    n2 = n_big;
                }                
            } else {
                n1 = n_big;
                n2 = n_small;
            }

            // std::cout<<"bounce: "<<j<<std::endl;
            // std::cout<<"intersected_big: "<<intersected_big<<", intersected_small: "<<intersected_small<<std::endl;
            // std::cout<<"in_big_circle: "<<in_big_circle<<std::endl;
            // std::cout<<"intersection_point: "<<intersection_point.x<<", "<<intersection_point.y<<std::endl;
            // std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;           

            //std::cout<<"cos_i "<<cos_i<<" n1: "<<n1<<", n2: "<<n2<<std::endl;
            double fresnel = fresnel_reflectance(cos_i, n1, n2);
            //double fresnel = 1.0;
            double random_number = random_double(0, 1);
            //std::cout<<"fresnel: "<<fresnel<<", random_number: "<<random_number<<std::endl;

            if (random_number < fresnel) {                
                // Reflect the ray
                reflect = true;
                // std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                // std::cout<<"normal: "<<normal.x<<", "<<normal.y<<std::endl;
                current_direction = current_direction - normal * (2 * current_direction.dot(normal));
                //std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                // // compute reflected angle and see if it equals to the incident angle
                // double reflected_angle = current_direction.dot(normal);
                // std::cout<<"incident angle "<<cos_i<<", reflected angle "<<reflected_angle<<std::endl;
                //abort();

                // fresnel = fresnel_reflectance(cos_i, n1, n2);
                // energy *= fresnel;

            } else {
                // Refract the ray
                reflect = false;
                // std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                // std::cout<<"normal: "<<normal.x<<", "<<normal.y<<std::endl;
                // std::cout<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;
                Vec2 refracted_direction = refract(current_direction, normal, n1, n2);
                // std::cout<<"refracted_direction: "<<refracted_direction.x<<", "<<refracted_direction.y<<std::endl;
                // abort();
                // Vec2 tmp_dir;
                // tmp_dir.x = -0.736996, tmp_dir.y = -0.675898;
                // Vec2 tmp_normal;
                // tmp_normal.x = 0.736996, tmp_normal.y = 0.675898;
                // Vec2 refracted_direction = refract(tmp_dir, tmp_normal, 1.0, 1.53);
                // if (std::isnan(refracted_direction.x) || std::isnan(refracted_direction.y)) {
                //     std::cout<<"current_direction: "<<current_direction.x<<", "<<current_direction.y<<std::endl;
                //     std::cout<<"normal: "<<normal.x<<", "<<normal.y<<std::endl;
                //     std::cout<<"n1 = "<<n1<<", n2 = "<<n2<<std::endl;
                //     abort();
                // }
                if (refracted_direction.norm() == 0) {
                    // Total internal reflection
                    current_direction = current_direction - normal * (2 * current_direction.dot(normal));
                } else {
                    current_direction = refracted_direction;
                }


                if (medium ==0 ){
                    medium = 1;
                } else if (medium == 1){
                    if (intersected_big){
                        medium = 0;
                    } else {
                        medium = 2;
                    }
                } else {                    
                    medium = 1;
                }
            }

            // Update the ray's origin to the intersection point
            previous_origin = current_origin;
            current_origin = intersection_point + current_direction * epsilon;

        }

        
        // Store the outgoing angle and energy if the ray interacted with the circle
        // std::cout<<"has_intersected: "<<has_intersected<<", exit: "<<exit<<std::endl;
        if (has_intersected && exit) {
            double outgoing_angle = std::atan2(current_direction.y, current_direction.x);
            #pragma omp critical
            {
                results.emplace_back(outgoing_angle, energy, y_start, current_origin, double(j));
                energies.push_back(energy);
                outgoing_angles.push_back(outgoing_angle);
            }
        }

        // // Store the outgoing angle and energy if the ray interacted with the circle
        // // std::cout<<"has_intersected: "<<has_intersected<<", exit: "<<exit<<std::endl;
        // if (has_intersected && exit && j==3) {
        //     double outgoing_angle = std::atan2(current_direction.y, current_direction.x);
        //     #pragma omp critical
        //     {
        //         results.emplace_back(outgoing_angle, energy, y_start, previous_origin, current_origin, double(j));
        //         energies.push_back(energy);
        //         outgoing_angles.push_back(outgoing_angle);
        //     }
        // }
    }

    // Save the energy distribution to a binary file    
    compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy, num_rays);
    // Save the accumulated energy to a binary file
    std::ofstream outfile(output_file, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << output_file << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<char*>(accumulated_energy.data()), num_bins * sizeof(double));
    outfile.close();
}



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <simulation_type> (single/multiple) [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --num_paths <int>         Number of paths (default: 1000000)" << std::endl;
        std::cerr << "  --radius_small <double>   Radius of small circles (default: 1.0)" << std::endl;
        std::cerr << "  --num_circles <int>       Number of small circles (default: 100)" << std::endl;
        std::cerr << "  --n2 <double>             Refractive index of the medium (default: 1.50)" << std::endl;
        std::cerr << "  --n3 <double>             Refractive index of the inner medium (default: 1.50)" << std::endl;
        std::cerr << "  --sigma_a <double>        Absorption coefficient (default: 0.0)" << std::endl;
        std::cerr << "  --max_bounces <int>       Maximum number of bounces (default: 50)" << std::endl;
        std::cerr << "  --num_bins <int>          Number of bins for angular distribution (default: 100)" << std::endl;
        std::cerr << "  --num_exp <int>           Number of experiments (default: 500)" << std::endl;
        std::cerr << "  --h <double>              h value (default: 0.4)" << std::endl;
        return 1;
    }

    std::string simulation_type = argv[1];

    double n1 = 1.0;
    unsigned int seed = 0;
    double radius_big = 10.0;

    // Default values
    int num_paths = 1000000;
    //int num_paths = 1;        
    double radius_small = 1.0;
    int num_circles = 100;
    double n2 = 1.50;
    double n3 = 1.50;
    double sigma_a = 0.0;
    int max_bounces = 100;
    int num_bins = 100;
    int num_exp = 500;
    double h = 0.4; // from -1 to 1

    std::cout<<"before parsing"<<std::endl;

    // Parse command-line arguments
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--num_paths" && i + 1 < argc) {
            num_paths = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--radius_small" && i + 1 < argc) {
            radius_small = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--num_circles" && i + 1 < argc) {
            num_circles = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--n2" && i + 1 < argc) {
            n2 = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--n3" && i + 1 < argc) {
            n3 = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--sigma_a" && i + 1 < argc) {
            sigma_a = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--max_bounces" && i + 1 < argc) {
            max_bounces = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--num_bins" && i + 1 < argc) {
            num_bins = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--num_exp" && i + 1 < argc) {
            num_exp = std::atoi(argv[i + 1]);
        } else if(std::string(argv[i]) == "--h" && i + 1 < argc) {
            h = std::atof(argv[i + 1]);
        }
    }
    

    // Debug output to confirm values
    std::cout << "Simulation type: " << simulation_type << std::endl;
    std::cout << "Number of paths: " << num_paths << std::endl;
    std::cout << "Radius of small circles: " << radius_small << std::endl;
    std::cout << "Number of small circles: " << num_circles << std::endl;
    std::cout << "Refractive index n2: " << n2 << std::endl;
    std::cout << "Refractive index n3: " << n3 << std::endl;
    std::cout << "Absorption coefficient sigma_a: " << sigma_a << std::endl;
    std::cout << "Maximum bounces: " << max_bounces << std::endl;
    std::cout << "Number of bins: " << num_bins << std::endl;  
    std::cout << "Number of experiments: " << num_exp << std::endl;
    std::cout << "h: " << h << std::endl;

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    if (simulation_type == "single") {
        //std::string output_file = "single.bin";
        // Create a stringstream to format the filename
        std::stringstream ss;
        ss << "oblique_single/single_sigma_" << std::fixed << std::setprecision(2) << sigma_a 
        << "_n2_" << std::fixed << std::setprecision(2) << n2 << ".bin";
        std::string output_file = ss.str();
        std::cout<<"output_file: "<<output_file<<std::endl;
        // Simulate rays intersecting with the large circle
        simulate_large_circle(num_paths, radius_big, n1, n2, sigma_a, seed, max_bounces, num_bins, output_file);
        //std::cout << "Large circle simulation completed.\n";
    } else if (simulation_type == "multiple") {
        //std::string output_file = "multiple_avg.bin";
        // Create a stringstream to format the filename
        std::stringstream ss;
        ss << "multiple_avg_sigma_" << std::fixed << std::setprecision(2) << sigma_a
        << "_n2_" << std::fixed << std::setprecision(2) << n2 << "_radius_"<<
        std::fixed << std::setprecision(1) << radius_small << ".bin";
        std::string output_file = ss.str();
        std::cout<<"output_file: "<<output_file<<std::endl;
        // Simulate rays interacting with small circles inside the large circle
        simulate_small_circles(num_paths, radius_big, radius_small, num_circles, n1, n2, sigma_a, seed, num_bins, num_exp, output_file);
        //std::cout << "Small circles simulation completed.\n";
     } else if (simulation_type == "multiple_fixedh") {
        //std::string output_file = "multiple_avg.bin";
        // Create a stringstream to format the filename
        std::stringstream ss;
        ss << "multiple_avg_sigma_" << std::fixed << std::setprecision(2) << sigma_a
        << "_n2_" << std::fixed << std::setprecision(2) << n2 << "_radius_"<<
        std::fixed << std::setprecision(1) << radius_small << "_h_"
        << std::fixed << std::setprecision(2) << h << ".bin";
        std::string output_file = ss.str();
        std::cout<<"output_file: "<<output_file<<std::endl;
        // Simulate rays interacting with small circles inside the large circle
        simulate_small_circles_h(num_paths, radius_big, radius_small, num_circles, n1, n2, sigma_a, seed, num_bins, num_exp, output_file, h);
        //std::cout << "Small circles simulation completed.\n";
    } else if (simulation_type == "double") {
        double sigma_a_big = 0.0;
        double sigma_a_small = 0.0;
        double n_big = n2;
        double n_small = n3;
        radius_small = 9.0;
        // Create a stringstream to format the filename
        std::stringstream ss;
        ss << "double_circles/multiple_avg_sigma_big_" << std::fixed << std::setprecision(2) << sigma_a_big
        << "_sigma_small_" << std::fixed << std::setprecision(2) << sigma_a_small
        << "_n_big_" << std::fixed << std::setprecision(2) << n_big 
        << "_n_small_" << std::fixed << std::setprecision(2) << n_small
        << "_radius_"<< std::fixed << std::setprecision(1) << radius_small << ".bin";
        std::string output_file = ss.str();
        std::cout<<"output_file: "<<output_file<<std::endl;
        // Simulate rays interacting with small circles inside the large circle

        simulate_double_circles(num_paths, radius_big, radius_small, n1, n_big, n_small, 
                                sigma_a_big, sigma_a_small, seed, max_bounces, num_bins, output_file);
    } else {
        std::cerr << "Invalid simulation type. Use 'single' or 'multiple'." << std::endl;
        return 1;
    }

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double> duration = end - start;

    // Output the elapsed time in seconds
    std::cout << "Simulation took " << duration.count() << " seconds." << std::endl;


    return 0;
}
