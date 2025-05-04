#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <random>

const bool debugMode = false;

const double G = 6.67430e-11; 

struct Vector3 {
    double x, y, z;

    Vector3 operator+(const Vector3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vector3 operator-(const Vector3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Vector3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vector3 normalize() const {
        double mag = magnitude();
        return {x / mag, y / mag, z / mag};
    }
};

struct Body {
    std::string name;
    double mass;
    Vector3 position;
    Vector3 velocity;
};

const std::unordered_map<std::string, double> radii = {
    {"Nap", 6.9634e8}, {"Merkúr", 2.4397e6}, {"Vénusz", 6.0518e6},
    {"Föld", 6.371e6}, {"Mars", 3.3895e6}, {"Jupiter", 6.9911e7},
    {"Szaturnusz", 5.8232e7}, {"Uránusz", 2.5362e7}, {"Neptunusz", 2.4622e7},
    {"Plútó", 1.1883e6}
};

const std::unordered_map<std::string, double> masses = {
    {"Nap", 1.989e30}, {"Merkúr", 3.301e23}, {"Vénusz", 4.867e24},
    {"Föld", 5.972e24}, {"Mars", 6.417e23}, {"Jupiter", 1.898e27},
    {"Szaturnusz", 5.683e26}, {"Uránusz", 8.681e25}, {"Neptunusz", 1.024e26},
    {"Plútó", 1.309e22}
};

double getRadius(const std::string& name) {
    auto it = radii.find(name);
    return (it != radii.end()) ? it->second : 1e6;
}

double getMass(const std::string& name) {
    auto it = masses.find(name);
    return (it != masses.end()) ? it->second : 1e20;
}

void computeForces(const std::vector<Body>& bodies, std::vector<Vector3>& forces) {
    size_t n = bodies.size();
    forces.assign(n, {0, 0, 0});

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            Vector3 diff = bodies[j].position - bodies[i].position;
            double dist = diff.magnitude();
            double forceMagnitude = G * bodies[i].mass * bodies[j].mass / (dist * dist);

            Vector3 force = diff.normalize() * forceMagnitude;
            forces[i] = forces[i] + force;
            forces[j] = forces[j] - force;

            if (debugMode) {
                std::cout << "Az erő " << bodies[i].name << " és " << bodies[j].name
                          << " között: " << forceMagnitude << " N, ennyi táv van köztük: "
                          << dist << " m\n";
            }
        }
    }
}

void rungeKuttaStep(std::vector<Body>& bodies, double dt) {
    size_t n = bodies.size();

    static std::vector<Vector3> k1_v(n), k1_p(n);
    static std::vector<Vector3> k2_v(n), k2_p(n);
    static std::vector<Vector3> k3_v(n), k3_p(n);
    static std::vector<Vector3> k4_v(n), k4_p(n);
    static std::vector<Vector3> forces(n);
    static std::vector<Body> tempBodies = bodies;

    computeForces(bodies, forces);
    for (size_t i = 0; i < n; ++i) {
        k1_v[i] = forces[i] * (dt / bodies[i].mass);
        k1_p[i] = bodies[i].velocity * dt;
    }

    for (size_t i = 0; i < n; ++i) {
        tempBodies[i].position = bodies[i].position + k1_p[i] * 0.5;
        tempBodies[i].velocity = bodies[i].velocity + k1_v[i] * 0.5;
    }
    computeForces(tempBodies, forces);
    for (size_t i = 0; i < n; ++i) {
        k2_v[i] = forces[i] * (dt / tempBodies[i].mass);
        k2_p[i] = tempBodies[i].velocity * dt;
    }

    for (size_t i = 0; i < n; ++i) {
        tempBodies[i].position = bodies[i].position + k2_p[i] * 0.5;
        tempBodies[i].velocity = bodies[i].velocity + k2_v[i] * 0.5;
    }
    computeForces(tempBodies, forces);
    for (size_t i = 0; i < n; ++i) {
        k3_v[i] = forces[i] * (dt / tempBodies[i].mass);
        k3_p[i] = tempBodies[i].velocity * dt;
    }

    for (size_t i = 0; i < n; ++i) {
        tempBodies[i].position = bodies[i].position + k3_p[i];
        tempBodies[i].velocity = bodies[i].velocity + k3_v[i];
    }
    computeForces(tempBodies, forces);
    for (size_t i = 0; i < n; ++i) {
        k4_v[i] = forces[i] * (dt / tempBodies[i].mass);
        k4_p[i] = tempBodies[i].velocity * dt;
    }

    for (size_t i = 0; i < n; ++i) {
        bodies[i].position = bodies[i].position + 
            (k1_p[i] + k2_p[i] * 2 + k3_p[i] * 2 + k4_p[i]) * (1.0 / 6.0);
        bodies[i].velocity = bodies[i].velocity + 
            (k1_v[i] + k2_v[i] * 2 + k3_v[i] * 2 + k4_v[i]) * (1.0 / 6.0);
    }
}

void initializeSolarSystem(std::vector<Body>& bodies, bool includeJupiter, unsigned int seed, const Vector3& externalVelocity) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> angleDist(0, 2 * M_PI);

    const std::unordered_map<std::string, double> orbitalRadii = {
        {"Merkúr", 5.79e10}, {"Vénusz", 1.082e11}, {"Föld", 1.496e11},
        {"Mars", 2.279e11}, {"Szaturnusz", 1.433e12},
        {"Uránusz", 2.877e12}, {"Neptunusz", 4.503e12}, {"Plútó", 5.906e12}
    };

    const std::unordered_map<std::string, double> orbitalSpeeds = {
        {"Merkúr", 47400}, {"Vénusz", 35020}, {"Föld", 29780},
        {"Mars", 24070}, {"Szaturnusz", 9680},
        {"Uránusz", 6800}, {"Neptunusz", 5430}, {"Plútó", 4740}
    };

    bodies = {
        {"Nap", getMass("Nap"), {0, 0, 0}, {0, 0, 0}}
    };

    for (const auto& [name, radius] : orbitalRadii) {
        double angle = angleDist(gen);
        double speed = orbitalSpeeds.at(name);

        Vector3 position = {radius * std::cos(angle), radius * std::sin(angle), 0};
        Vector3 velocity = {-speed * std::sin(angle), speed * std::cos(angle), 0};

        bodies.push_back({name, getMass(name), position, velocity});
    }

    if (includeJupiter) {
        double angle = angleDist(gen);
        double radius = 7.785e11;
        double speed = 13070;

        Vector3 position = {radius * std::cos(angle), radius * std::sin(angle), 0};
        Vector3 velocity = {-speed * std::sin(angle), speed * std::cos(angle), 0};

        bodies.push_back({"Jupiter", getMass("Jupiter"), position, velocity});
    }

    bodies.push_back({"IdegenObjektum", 1e12, {6e12, 0, 0}, externalVelocity});
}

bool detectCollision(const std::vector<Body>& bodies, int& earthCollisionCount) {
    const double externalRadius = 1e6;

    for (const auto& body : bodies) {
        if (body.name == "IdegenObjektum") {
            for (const auto& other : bodies) {
                if (other.name != "IdegenObjektum") {
                    const double otherRadius = getRadius(other.name);
                    const double collisionDistance = (externalRadius + otherRadius) * 1.5;

                    Vector3 diff = body.position - other.position;
                    double distance = diff.magnitude();

                    if (distance < collisionDistance) {
                        if (other.name == "Föld") {
                            earthCollisionCount++;
                            if (debugMode) {
                                std::cout << "Ütközés az IdegenObjektum és a Föld között! Hibás elrendezés!\n"
                                          << "Táv: " << distance << " m\n";
                            }
                        } else if (debugMode) {
                            std::cout << "Ütközés az IdegenObjektum és a(z) " << other.name
                                      << " között! Hibás elrendezés! Táv: " << distance << " m\n";
                        }
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void runSimulation(bool includeJupiter, int numSimulations, double dt, int steps, unsigned int seed, int& totalCollisions) {
    totalCollisions = 0;

    for (int sim = 1; sim <= numSimulations; ++sim) {
        std::vector<Body> bodies;

        int velocityIterations = debugMode ? 1 : 1000;

        for (int i = 0; i < velocityIterations; ++i) {
            double angle = 2 * M_PI * i / velocityIterations;
            Vector3 externalVelocity = {72000 * std::cos(angle), 72000 * std::sin(angle), 0};

            initializeSolarSystem(bodies, includeJupiter, seed, externalVelocity);

            if (debugMode) {
                std::cout << "Szimuláció " << sim << ", sebesség szöge: " << angle << " rad\n";
                for (const auto& body : bodies) {
                    std::cout << "Objektum: " << body.name
                              << ", Pozíció: (" << body.position.x << ", " << body.position.y << ", " << body.position.z << ")"
                              << ", Sebesség: (" << body.velocity.x << ", " << body.velocity.y << ", " << body.velocity.z << ")\n";
                }
            }

            for (int step = 0; step < steps; ++step) {
                rungeKuttaStep(bodies, dt);

                if (debugMode) {
                    std::cout << "A " << step << " lépéssorszám kész.\n";
                    for (const auto& body : bodies) {
                        std::cout << "Objektum: " << body.name
                                  << ", Pozíció: (" << body.position.x << ", " << body.position.y << ", " << body.position.z << ")"
                                  << ", Sebesség: (" << body.velocity.x << ", " << body.velocity.y << ", " << body.velocity.z << ")\n";
                    }
                }

                if (detectCollision(bodies, totalCollisions)) {
                    if (debugMode) {
                        std::cout << "A szimuláció megállt hibás elrendezés miatt a  " << step << " lépésben.\n";
                    }
                    break;
                }
            }

            if (debugMode) {
                break;
            }
        }

        if (debugMode) {
            break;
        }
    }
}

int main() {
    int numSimulations = debugMode ? 1 : 1;
    double dt = 3600; 
    int steps = debugMode ? 1 : 24 * 365 * 5; 
    unsigned int seed = 42; 

    int collisionsWithJupiter = 0;
    int collisionsWithoutJupiter = 0;

    std::cout << "\n\n\nJupiteres szimuláció...\n\n\n";
    runSimulation(true, numSimulations, dt, steps, seed, collisionsWithJupiter);

    std::cout << "\n\n\nJupiter nélküli szimuláció...\n\n\n";
    runSimulation(false, numSimulations, dt, steps, seed, collisionsWithoutJupiter);

    if (!debugMode) {
        std::cout << "Ütközések száma Jupiterrel: " << collisionsWithJupiter << "\n";
        std::cout << "Ütközések száma Jupiter nélkül: " << collisionsWithoutJupiter << "\n";
    }

    return 0;
}