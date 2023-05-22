#include <SFML/Graphics.hpp>
#include <vector>
#include <thread>
#include <random>
#include <cmath>
#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Simple struct to hold particle data
struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    float orientation;
    float angularVelocity;
    sf::Vector2f force;
    float torque;
    int type;
    //std::vector<Particle*> linked;
    //Having view relating object in the model object is not ideal
    //but in SFML not rebuilding the graphical object is significantly faster
    sf::CircleShape shape;
};

sf::Vector2f multiply(const sf::Vector2f& vec, float scalar) {
    return sf::Vector2f(vec.x * scalar, vec.y * scalar);
}

float angleFromVector(sf::Vector2f v) {
    return std::atan2(v.y, v.x);
}

static int K_PARTICLES = 100;
static int WORLD_WIDTH = 1920;
static int WORLD_HEIGTH = 1080;
static int DOT_SIZE = 2;
static sf::Vector2f DOT_OFSET = sf::Vector2f(DOT_SIZE, DOT_SIZE);

float norm(const sf::Vector2f& vec) {
    return std::sqrt(vec.x * vec.x + vec.y * vec.y);
}

void normalize(sf::Vector2f& vec)
{
    float normVec = norm(vec);
    if (!normVec == 0) vec /= normVec;
}

sf::Color getColor(int type) {
    // Convert the type to a hue value between 0 and 360 degrees
    float hue = (type % 16) * (360.0f / 16.0f);

    // For simplicity, we'll keep saturation and lightness constant
    float saturation = 1.0f;
    float lightness = 0.5f;

    // Convert from HSL to RGB using the formula
    float c = (1.0f - std::abs(2.0f * lightness - 1.0f)) * saturation;
    float x = c * (1.0f - std::abs(std::fmod(hue / 60.0f, 2.0f) - 1.0f));
    float m = lightness - c / 2.0f;

    float r = 0, g = 0, b = 0;
    if (0 <= hue && hue < 60) {
        r = c, g = x, b = 0;
    } else if (60 <= hue && hue < 120) {
        r = x, g = c, b = 0;
    } else if (120 <= hue && hue < 180) {
        r = 0, g = c, b = x;
    } else if (180 <= hue && hue < 240) {
        r = 0, g = x, b = c;
    } else if (240 <= hue && hue < 300) {
        r = x, g = 0, b = c;
    } else if (300 <= hue && hue < 360) {
        r = c, g = 0, b = x;
    }

    return sf::Color((r + m) * 255, (g + m) * 255, (b + m) * 255);
}

struct Model {
    std::vector<Particle> particles;
    std::random_device rd;
    std::mt19937 gen;

    Model() : gen(rd())
    {
        init();
    }

    void init()
    {
        // Initialize random number generator
        std::uniform_int_distribution<> disType(0, 15);
        std::uniform_real_distribution<> disX(-WORLD_WIDTH/16, WORLD_WIDTH/16);
        std::uniform_real_distribution<> disY(-WORLD_HEIGTH/16, WORLD_HEIGTH/16);
        std::uniform_real_distribution<> disV(-2, 2);
        std::uniform_real_distribution<> disA(0, 2.0*M_PI);

        // Initialize particles
        particles.clear();
        for (int i = 0; i < K_PARTICLES; ++i) {
            Particle p;
            p.position = sf::Vector2f(disX(gen), disY(gen));
            p.velocity = sf::Vector2f(disV(gen), disV(gen));
            p.orientation = disA(gen);
            p.angularVelocity = 0.0;
            p.type = disType(gen);
            p.shape = sf::CircleShape(DOT_SIZE);
            p.shape.setOrigin(DOT_SIZE, DOT_SIZE);
            p.shape.setPosition(p.position);
            particles.push_back(p);
        }
    }

    void calculateForceAndTorque(const sf::Vector2f& r, float rNorm, float orientation1, float orientation2, sf::Vector2f& force, float& torque)
    {
        // The interaction strength
        float strengthF = 10.0f;
        float strengthT = 1.0f;

        // Calculate the relative orientation
        float deltaOrientation = orientation2 - orientation1;
        float cosOrientation = std::cos(deltaOrientation);
        float sinOrientation = std::sin(deltaOrientation);
        float dirAngle = angleFromVector(r);
        float deltaForceOrientation = orientation1-dirAngle;
        float sinDeltaForceOrientation = std::sin(deltaForceOrientation);
        // Make sure deltaOrientation is between -pi and pi
        //while (deltaOrientation > M_PI) deltaOrientation -= 2 * M_PI;
        //while (deltaOrientation < -M_PI) deltaOrientation += 2 * M_PI;

        // Calculate the force
        // The force magnitude depends on the relative orientation
        float forceMagnitude = strengthF * cosOrientation * (sinDeltaForceOrientation * sinDeltaForceOrientation - 0.2);
        // Solid repulsion
        if (rNorm < 2.0*DOT_SIZE)
            forceMagnitude = -100.0;

        force = r * (forceMagnitude / rNorm);

        // Calculate the torque
        // The torque magnitude depends on the relative orientation and the distance
        //torque = strengthT * sinOrientation * rNorm;
    }

    void step()
    {
        //const float attractionStrength = 0.05f;
        const float interactionRadius = 30.0f;
        //const float attractionMinThreshold = 5.0f;
        //const float linkingThreshold = 30.0f;

        // Calculate the force and torque on particle p due to all other particles
        for (auto& p : particles) {
            p.force = sf::Vector2f(0.0, 0.0);
            p.torque = 0.0;
            for (Particle& other : particles) {
                if (&other != &p) {  // Avoid self-interaction
                    sf::Vector2f r = other.position - p.position;
                    float rNorm = norm(r);
                    if (rNorm < interactionRadius) {  // Consider only particles within the interaction radius
                        // Calculate force and torque using appropriate model
                        sf::Vector2f force(0.0, 0.0);
                        float torque = 0.0;
                        calculateForceAndTorque(r, rNorm, p.orientation, other.orientation, force, torque);
                        p.force += force;
                        p.torque += torque;
                    }
                }
                // Floc behavior
                //p.velocity = multiply(p.velocity, 0.999);
                //p.velocity += multiply(other.velocity, 0.0001);
                //p.angularVelocity = p.angularVelocity * 0.999;
                //p.angularVelocity += other.angularVelocity * 0.0001;
            }

            // Containing forces
            if (p.position.x > WORLD_WIDTH/2) p.force += sf::Vector2f(-0.1,0.0);
            if (p.position.x < -WORLD_WIDTH/2) p.force += sf::Vector2f(0.1,0.0);
            if (p.position.y > WORLD_HEIGTH/2) p.force += sf::Vector2f(0.0,-0.1);
            if (p.position.y < -WORLD_HEIGTH/2) p.force += sf::Vector2f(0.0,0.1);
        }

        for (auto& p : particles) {
            // Calculate the acceleration and angular acceleration (assuming mass and moment of inertia = 1)
            sf::Vector2f acceleration = p.force;
            float angularAcceleration = p.torque;
            float dt = 0.01;

            // Using naive algo
            p.velocity = p.velocity + acceleration * dt;
            p.angularVelocity = p.angularVelocity + angularAcceleration * dt;
            float maxVelocity = 50.0;
            if (norm(p.velocity) > maxVelocity) {
                normalize(p.velocity);
                p.velocity *= maxVelocity;
            }
            p.velocity -= multiply(p.position, 0.001);
            p.velocity = multiply(p.velocity, 0.99);
            float maxAngularVelocity = 10.0;
            if (p.angularVelocity > maxAngularVelocity) p.angularVelocity = maxAngularVelocity;
            if (p.angularVelocity < -maxAngularVelocity) p.angularVelocity = -maxAngularVelocity;
            p.position = p.position + p.velocity * dt;
            p.orientation = p.orientation + p.angularVelocity * dt;

            // Update position and velocity using Velocity Verlet algorithm
            //p.position = p.position + p.velocity * dt + 0.5f * acceleration * dt * dt;
            //sf::Vector2f newAcceleration = p.force;  // Recalculate force here if necessary
            //p.velocity = p.velocity + 0.5f * (acceleration + newAcceleration) * dt;

            // Update orientation and angular velocity
            //p.orientation = p.orientation + p.angularVelocity * dt + 0.5f * angularAcceleration * dt * dt;
            //float newAngularAcceleration = p.torque;  // Recalculate torque here if necessary
            //p.angularVelocity = p.angularVelocity + 0.5f * (angularAcceleration + newAngularAcceleration) * dt;

            // Update the particle's shape position
            p.shape.setPosition(p.position);
            p.shape.setFillColor(getColor(p.type));
        }

            // sort by distance
            //std::sort(neighbors.begin(), neighbors.end());

            //int nCount = 0;
            //int dCut = 8;
            //for (auto& pair : neighbors) {
            //    nCount++;
            //    float distance = pair.first;
            //    Particle* other = pair.second;

            //    sf::Vector2f delta = other->position - p.position;
            //    delta /= distance; // Normalize delta
            //    float f = 1.0*(dCut-(nCount-1.0))/dCut;
            //    if (f < -0.2) f = -0.2;
            //    if (distance > attractionMinThreshold) {
            //        p.velocity += f * delta * attractionStrength;
            //    } else {
            //        float rep = 0.2;
            //        p.velocity -= rep * delta * attractionStrength;
            //    }
                //if (nCount >= dCut){
                //    p.velocity -= delta * attractionStrength;
                //}

                // Floc behavior
            //    p.velocity = multiply(p.velocity, 0.999);
            //    p.velocity += multiply(other->velocity, 0.0001);
            //}

            //std::uniform_real_distribution<> disPV(-0.1, 0.1);
            //p.velocity += sf::Vector2f(disPV(gen),disPV(gen));

            //p.velocity -= multiply(p.position, 0.001);

            //if (norm(p.velocity) < 1.0) {
            //    normalize(p.velocity);
            //}

            //if (norm(p.velocity) > 2.0) {
            //    p.velocity = multiply(p.velocity, 0.9);
            //}

            // Apply force
            //for (auto& other : p.linked) {
            //    sf::Vector2f delta = other->position - p.position;
            //    float distance = std::sqrt(delta.x*delta.x + delta.y*delta.y);
            //    delta /= distance; // Normalize delta
            //    p.velocity += delta * attractionStrength;
            //}
            //normalize(p.velocity);

            // TODO Torus world is problematic to draw lines between points
            // Torus world
            //if (p.position.x > WORLD_WIDTH/2) p.position.x -= WORLD_WIDTH;
            //if (p.position.x < -WORLD_WIDTH/2) p.position.x += WORLD_WIDTH;
            //if (p.position.y > WORLD_HEIGTH/2) p.position.y -= WORLD_HEIGTH;
            //if (p.position.y < -WORLD_HEIGTH/2) p.position.y += WORLD_HEIGTH;
    }
};

void drawModel(sf::RenderWindow& ioWindow, const Model& iModel) {
    for (const auto& p : iModel.particles)
    {
        /*for (auto& other : p.linked) {
            sf::VertexArray lines(sf::Lines, 2);
            lines[0].position = p.position + DOT_OFSET;
            lines[1].position = other->position + DOT_OFSET;
            lines[0].color = sf::Color::White;
            lines[1].color = sf::Color::White;
            ioWindow.draw(lines);
        }*/

        float tailLength = 10.0;
        sf::Vertex line[] =
        {
            sf::Vertex(p.position, sf::Color::White),
            sf::Vertex(p.position + tailLength*sf::Vector2f(std::cos(p.orientation), std::sin(p.orientation)), sf::Color::White)
        };

        //Force debug
        //sf::Vertex lineF[] =
        //{
        //    sf::Vertex(p.position, sf::Color::Red),
        //    sf::Vertex(p.position + p.force, sf::Color::Red)
        //};
        //ioWindow.draw(lineF, 2, sf::Lines);

        ioWindow.draw(line, 2, sf::Lines);
        ioWindow.draw(p.shape);
    }
}

// Main function
int main()
{
    // Test
    //std::cout << angleFromVector(sf::Vector2f(-1.0,-2.0)) << std::endl;

    // Create the main window
    sf::RenderWindow window(sf::VideoMode(WORLD_WIDTH, WORLD_HEIGTH), "Particle system");

    // Model
    Model myModel;

    // Create a view with the same size as the window
    sf::View view(sf::FloatRect(-WORLD_WIDTH/2, -WORLD_HEIGTH/2, WORLD_WIDTH, WORLD_HEIGTH));

    // Variables to store the state of mouse dragging
    bool isDragging = false;
    sf::Vector2f lastMousePos;

    // main loop
    while (window.isOpen()) {
        // Handle events
        sf::Event event;
        while (window.pollEvent(event)) {
            switch (event.type)
            {
            case sf::Event::Closed:
                window.close();
                break;
            case sf::Event::MouseWheelScrolled:
                // If the wheel scrolled up, zoom in, else zoom out
                if (event.mouseWheelScroll.delta > 0)
                    view.zoom(0.9f);
                else
                    view.zoom(1.1f);
                break;
            case sf::Event::MouseButtonPressed:
                if (event.mouseButton.button == sf::Mouse::Left) {
                    isDragging = true;
                    lastMousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                }
                break;
            case sf::Event::MouseButtonReleased:
                if (event.mouseButton.button == sf::Mouse::Left)
                    isDragging = false;
                break;
            }
        }

        // Check for mouse dragging
        if (isDragging && sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            const sf::Vector2f mousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));
            const sf::Vector2f delta = lastMousePos - mousePos;

            view.move(delta);
        }

        // Refresh
        window.setView(view);

        // Update the last mouse position
        lastMousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));

        // Model update
        myModel.step();

        // Clear screen
        window.clear();

        // Draw particles
        drawModel(window, myModel);

        // Update the window
        window.display();
    }

    return 0;
}
