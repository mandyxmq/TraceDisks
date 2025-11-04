#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <sstream>

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

double random_double(double min, double max) {
    static thread_local std::mt19937 gen_local;
    static thread_local bool seeded = false;
    if (!seeded) {
        unsigned int thread_id = omp_get_thread_num();
        gen_local.seed(1234 + thread_id); // base seed + thread ID
        seeded = true;
    }
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen_local);
}

// Function to create random circles using grid-based optimization
std::vector<Circle> create_random_circles(double radius_big, double radius_small, int num_circles, double epsilon) {
    std::vector<Circle> circles;
    int attempts = 0;
    int max_attempts = 1000;

    // This is the correct cell size
    double cell_size = 2 * radius_small; 

    // The grid must cover the entire diameter (2 * radius_big)
    int grid_size = static_cast<int>((2 * radius_big) / cell_size) + 1;
    std::vector<std::vector<std::vector<Circle>>> grid(grid_size, std::vector<std::vector<Circle>>(grid_size));

    // The lambda should map coordinates from [-radius_big, +radius_big] 
    // to the grid index [0, grid_size - 1]
    auto get_grid_index = [&](double x, double y) {
        // (x + radius_big) shifts the range from [-R, +R] to [0, 2R]
        int ix = static_cast<int>((x + radius_big) / cell_size);
        int iy = static_cast<int>((y + radius_big) / cell_size);
        
        // Use clamp to be safe
        ix = std::clamp(ix, 0, grid_size - 1);
        iy = std::clamp(iy, 0, grid_size - 1);
        
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

bool intersect_ray_circle(const Vec2& ray_origin, const Vec2& ray_direction, 
    const Vec2& circle_center, double circle_radius, Vec2& intersection_point) {

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

void compute_distributions(const std::vector<double>& outgoing_angles, const std::vector<double>& energies, 
    int num_bins, std::vector<double>& accumulated) {

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

bool trace_ray(const Vec2& start_point, Vec2 ray_direction, 
            const std::vector<Circle>& circles, double radius_big, 
            double n1, double n2, double sigma_a, double incident_angle, 
            std::vector<double>& outgoing_angles, std::vector<double>& energies) {

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
        double transmission = std::exp(-sigma_a * entry_to_intersection_length / (2 * radius_big));
        energy *= transmission;

        if (energy < lowest_energy) {
            return interacted; // Path terminated with negligible energy
        }

        double cos_i = -ray_direction.dot(closest_normal);
        double fresnel = fresnel_reflectance(cos_i, n1, n2);
        double random_number = random_double(0, 1);

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
    }

    if (interacted) {
        double outgoing_angle = std::atan2(ray_direction.y, ray_direction.x) - incident_angle;
        outgoing_angle = fmod(outgoing_angle + 2 * M_PI, 2 * M_PI);
        outgoing_angles.push_back(outgoing_angle);
        energies.push_back(energy);
    }

    return interacted;
}

bool trace_ray_record(const Vec2& start_point, Vec2 ray_direction,
                    const std::vector<Circle>& circles, double radius_big,
                    double n1, double n2, double sigma_a, double incident_angle,
                    std::vector<double>& outgoing_angles, std::vector<double>& energies,
                    std::vector<Vec2>* out_points // may be nullptr
) {
    double energy = 1.0;
    Vec2 ray_origin = start_point;
    bool interacted = false;
    double lowest_energy = 1e-4;
    double epsilon = 1e-6;

    if (out_points) out_points->push_back(ray_origin);

    for (int i = 0; i < 100; ++i) {
        Vec2 closest_intersection;
        Circle closest_circle = {{0,0}, 0};
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

        if (!intersection_found) break;

        interacted = true;
        double entry_to_intersection_length = (closest_intersection - ray_origin).norm();
        double transmission = std::exp(-sigma_a * entry_to_intersection_length / (2 * radius_big));
        energy *= transmission;

        if (energy < lowest_energy) return interacted;

        double cos_i = -ray_direction.dot(closest_normal);
        double fresnel = fresnel_reflectance(cos_i, n1, n2);
        double rnd = random_double(0, 1);

        if (out_points) out_points->push_back(closest_intersection);

        if (rnd < fresnel) {
            ray_direction = reflect(ray_direction, closest_normal);
        } else {
            Vec2 refracted_direction = refract(ray_direction, closest_normal, n1, n2);
            if (refracted_direction.norm() == 0) {
                ray_direction = reflect(ray_direction, closest_normal);
            } else {
                ray_direction = refracted_direction;
            }
        }
        ray_origin = closest_intersection;
    }

    if (interacted) {
        double outgoing_angle = std::atan2(ray_direction.y, ray_direction.x) - incident_angle;
        outgoing_angle = fmod(outgoing_angle + 2 * M_PI, 2 * M_PI);
        outgoing_angles.push_back(outgoing_angle);
        energies.push_back(energy);
    }
    return interacted;
}

void simulate_multiple_circles(int num_paths, double radius_big, double radius_small, 
                            int num_circles, double n1, double n2, double sigma_a, 
                            int num_bins, int num_exp, double epsilon, std::string output_file) {
    double total_energy = 0.0;
    double incident_angle = 0.0; 
    std::vector<double> accumulated_energy(num_bins, 0.0);
    
    for (int index = 0; index < num_exp; ++index){
        // Print progress as "done/total" (e.g. "42/500") and overwrite in place
        std::cout << "\rExperiment " << (index + 1) << "/" << num_exp << std::flush;
        std::vector<Circle> circles = create_random_circles(radius_big, radius_small, num_circles, epsilon);

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

            #pragma omp for
            for (int i = 0; i < num_paths; ++i) {
                double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_paths));  
                Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Start just left of the circle
                Vec2 ray_direction = {std::cos(incident_angle), std::sin(incident_angle)};
                
                double theta = (M_PI - 2 * M_PI / 100.0) / num_paths * i + M_PI / 2 + M_PI / 100.0;
                
                // Trace rays and store results in local vectors
                bool interacted = trace_ray_record(ray_origin, ray_direction, circles, radius_big, n1, n2, sigma_a, 
                              incident_angle, local_outgoing_angles, local_energies, nullptr);
               
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
            }
        }

        // Save the energy distribution to a binary file
        compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy);
    }

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
}

void simulate_multiple_circles_plot(int num_paths, double radius_big, double radius_small,
                                int num_circles, double n1, double n2, double sigma_a,
                                int max_traces, double epsilon,
                                const std::string& circles_csv, const std::string& rays_json) {

    auto circles = create_random_circles(radius_big, radius_small, num_circles, epsilon);

    // Write circles.csv: cx,cy,r (plus the big boundary as radius_big at 0,0)
    {
        std::ofstream fc(circles_csv);
        fc << std::setprecision(17);
        fc << "cx,cy,r\n";
        // big boundary (optional: first row)
        fc << 0.0 << "," << 0.0 << "," << radius_big << "\n";
        for (const auto& c : circles) {
            fc << c.center.x << "," << c.center.y << "," << c.radius << "\n";
        }
    }

    // A small sample of rays to record
    std::ofstream fr(rays_json);
    fr << std::setprecision(17);

    double incident_angle = 0.0; // horizontal rays from left
    int recorded = 0;
    int stride = std::max(1, num_paths / std::max(1, max_traces)); // spread them out

    // We’ll do this serially to keep the file output simple/deterministic.
    for (int i = 0; i < num_paths && recorded < max_traces; i += stride) {
        double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_paths));
        Vec2 ray_origin = {-radius_big - 1.0, y_start};
        Vec2 ray_direction = {std::cos(incident_angle), std::sin(incident_angle)};

        std::vector<double> angles, energies;
        std::vector<Vec2> points;
        bool interacted = trace_ray_record(
            ray_origin, ray_direction,
            circles, radius_big, n1, n2, sigma_a, incident_angle,
            angles, energies, &points
        );

        if (!interacted) continue;

        // Write one ray as a single JSON object per line (JSON)
        fr << "{";
        fr << "\"y_start\":" << y_start << ",\"points\":[";
        for (size_t k = 0; k < points.size(); ++k) {
            fr << "[" << points[k].x << "," << points[k].y << "]";
            if (k + 1 < points.size()) fr << ",";
        }
        fr << "]";
        if (!energies.empty())
            fr << ",\"energy\":" << energies.back();
        if (!angles.empty())
            fr << ",\"out_angle\":" << angles.back();
        fr << "}\n";

        recorded++;
    }
    fr.close();

    std::cout << "Plot data saved to " << circles_csv << " and " << rays_json << std::endl;
}

void simulate_single_circle(int num_rays, double radius_big, double n1, double n2, double sigma_a, int max_bounces, int num_bins, std::string output_file) {
    Vec2 circle_center = {0.0, 0.0};  // Center of the big circle
    std::vector<double> energies;
    std::vector<double> outgoing_angles;
    double min_energy = 1e-4;
    double epsilon = 1e-6;
    std::vector<double> accumulated_energy(num_bins, 0.0);

    #pragma omp parallel for
    for (int i = 0; i < num_rays; ++i) {
        // Uniformly distribute starting points
        double y_start = -radius_big + (2.0 * radius_big * double(i) / double(num_rays)) + 1.0 * radius_big / double(num_rays); 
        Vec2 ray_origin = {-radius_big - 1.0, y_start};  // Rays coming from the left of the circle
        Vec2 ray_direction = {1.0, 0.0};  // Horizontal ray moving rightwards
        
        Vec2 current_origin = ray_origin;
        Vec2 previous_origin = ray_origin;
        Vec2 current_direction = ray_direction;
        double energy = 1.0;  // Initial energy of each ray
        bool has_intersected = false;
        bool exit = false;
        int j = 0;

        for (j = 0; j < max_bounces; ++j) {  
            Vec2 intersection_point;
            if (!intersect_ray_circle(current_origin, current_direction, circle_center, radius_big, intersection_point)) {
                exit = true;
                break;  // No intersection with the circle; ray exits
            }

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
            double random_number = random_double(0, 1);
            if (random_number < fresnel) {
                // Reflect the ray
                current_direction = current_direction - normal * (2 * current_direction.dot(normal));
            } else {
                // Refract the ray
                Vec2 refracted_direction = refract(current_direction, normal, n1, n2);
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
        if (has_intersected && exit) {
            double outgoing_angle = std::atan2(current_direction.y, current_direction.x);
            #pragma omp critical
            {
                energies.push_back(energy);
                outgoing_angles.push_back(outgoing_angle);
            }
        }
    }

    // Save the energy distribution to a binary file    
    compute_distributions(outgoing_angles, energies, num_bins, accumulated_energy);
    // normalize accumulated_energy
    for (int i = 0; i < num_bins; ++i) {
        accumulated_energy[i] /= num_rays;
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <simulation_type> (single/multiple) [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --num_paths <int>         Number of paths (default: 1000000)" << std::endl;
        std::cerr << "  --radius_small <double>   Radius of small circles (default: 1.0)" << std::endl;
        std::cerr << "  --num_circles <int>       Maximum number of small circles (default: 100)" << std::endl;
        std::cerr << "  --n1 <double>             Refractive index of the environment (default: 1.00)" << std::endl;
        std::cerr << "  --n2 <double>             Refractive index of the small circles (default: 1.50)" << std::endl;
        std::cerr << "  --sigma_a <double>        Absorption coefficient (default: 0.0)" << std::endl;
        std::cerr << "  --max_bounces <int>       Maximum number of bounces (default: 50)" << std::endl;
        std::cerr << "  --num_bins <int>          Number of bins for angular distribution (default: 100)" << std::endl;
        std::cerr << "  --num_exp <int>           Number of experiments (default: 100)" << std::endl;
        std::cerr << "  --epsilon <double>        Smallest gap between small circles (default: 0.1)" << std::endl;
        return 1;
    }

    std::string simulation_type = argv[1];

    // Default values
    double n1 = 1.0;
    double radius_big = 10.0;
    int num_paths = 1000000;
    double radius_small = 1.0;
    int num_circles = 100;
    double n2 = 1.50;
    double n3 = 1.50;
    double sigma_a = 0.0;
    int max_bounces = 100;
    int num_bins = 100;
    int num_exp = 100;
    double epsilon = 0.1;
    bool plot = false;

    // Parse command-line arguments
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--num_paths" && i + 1 < argc) {
            num_paths = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--radius_small" && i + 1 < argc) {
            radius_small = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--num_circles" && i + 1 < argc) {
            num_circles = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--n1" && i + 1 < argc) {
            n1 = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--n2" && i + 1 < argc) {
            n2 = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--sigma_a" && i + 1 < argc) {
            sigma_a = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--max_bounces" && i + 1 < argc) {
            max_bounces = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--num_bins" && i + 1 < argc) {
            num_bins = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--num_exp" && i + 1 < argc) {
            num_exp = std::atoi(argv[i + 1]);
        } else if (std::string(argv[i]) == "--epsilon" && i + 1 < argc) {
            epsilon = std::atof(argv[i + 1]);
        } else if (std::string(argv[i]) == "--plot") {
            plot = true;
        }
    }
    

    // Debug output to confirm values
    std::cout << "Simulation type: " << simulation_type << std::endl;
    std::cout << "Number of paths: " << num_paths << std::endl;
    std::cout << "Radius of small circles: " << radius_small << std::endl;
    std::cout << "Maximum number of small circles: " << num_circles << std::endl;
    std::cout << "Refractive index outside: " << n1 << std::endl;
    std::cout << "Refractive index of the small circles: " << n2 << std::endl;
    std::cout << "Absorption coefficient sigma_a: " << sigma_a << std::endl;
    std::cout << "Maximum bounces: " << max_bounces << std::endl;
    std::cout << "Number of bins: " << num_bins << std::endl;  
    if (simulation_type == "multiple"){
        if (plot) {
            std::cout << "Path recording enabled." << std::endl;
        } else {
            std::cout << "Number of experiments: " << num_exp << std::endl;
        }
        std::cout << "Smallest gap between small circles: " << epsilon << std::endl;
    }

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    if (simulation_type == "single") {
        // Create a stringstream to format the filename
        std::stringstream ss;
        ss << "single_sigma_" << std::fixed << std::setprecision(2) << sigma_a 
        << "_n1_" << std::fixed << std::setprecision(2) << n1
        << "_n2_" << std::fixed << std::setprecision(2) << n2 << ".bin";
        std::string output_file = ss.str();
        std::cout<<"output_file: "<<output_file<<std::endl;
        // Simulate rays intersecting with the large circle
        simulate_single_circle(num_paths, radius_big, n1, n2, sigma_a, max_bounces, num_bins, output_file);

    } else if (simulation_type == "multiple") {
        if (plot) {
            simulate_multiple_circles_plot(
                num_paths, radius_big, radius_small,
                num_circles, n1, n2, sigma_a,
                /*max_traces=*/30, epsilon,
                "circles.csv", "rays.json"
            );
        } else {
            std::cout << "Running multiple circles simulation..." << std::endl;
            // Create a stringstream to format the filename
            std::stringstream ss;
            ss << "multiple_avg_sigma_" << std::fixed << std::setprecision(2) << sigma_a
            << "_n1_" << std::fixed << std::setprecision(2) << n1
            << "_n2_" << std::fixed << std::setprecision(2) << n2 << "_radius_"<<
            std::fixed << std::setprecision(1) << radius_small << ".bin";
            std::string output_file = ss.str();
            std::cout<<"output_file: "<<output_file<<std::endl;

            // Simulate rays interacting with small circles inside the large circle
            simulate_multiple_circles(num_paths, radius_big, radius_small, num_circles, n1, n2, sigma_a, num_bins, num_exp, epsilon, output_file);
        }

     } else {
        std::cerr << "Invalid simulation type. Use 'single' or 'multiple'." << std::endl;
        return 1;
    }

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double> duration = end - start;

    // Output the elapsed time in seconds
    std::cout << std::endl;
    std::cout << "Simulation took " << duration.count() << " seconds." << std::endl;

    return 0;
}
