#include <SFML/Graphics.hpp>
#include <vector>
#include <thread>
#include <random>
#include <cmath>
#include <iostream>
#include <cassert>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static int K_PARTICLES = 64;
static int WORLD_WIDTH = 1920;
static int WORLD_HEIGTH = 1080;
static int DOT_SIZE = 2;
static int K_NB_TYPE = 16;
static sf::Vector2f DOT_OFSET = sf::Vector2f(DOT_SIZE, DOT_SIZE);

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

sf::Vector2f unitVectorFromAngle(float iAngle) {
    return sf::Vector2f(std::cos(iAngle), std::sin(iAngle));
}

float colinearFactor(sf::Vector2f v1, sf::Vector2f v2) {
    float dotProduct = v1.x * v2.x + v1.y * v2.y;
    float lengthsProduct = std::sqrt(v1.x * v1.x + v1.y * v1.y) * std::sqrt(v2.x * v2.x + v2.y * v2.y);

    // prevent division by zero
    if(lengthsProduct == 0.f)
        return 0.f;

    float cosineOfAngle = dotProduct / lengthsProduct;

    // Clamp the value to the [-1, 1] range, in case of numerical instability
    cosineOfAngle = std::max(-1.f, std::min(1.f, cosineOfAngle));

    return cosineOfAngle;
}

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
    float hue = (type % K_NB_TYPE) * (360.0f / 16.0f);

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
    std::vector<Particle> particles; //BAD choice of container if dynamic
    std::random_device rd;
    std::mt19937 gen;

    Model() : gen(rd())
    {
        init();
    }

    void init()
    {
        // Initialize random number generator
        std::uniform_int_distribution<> disType(0, K_NB_TYPE-1);
        std::uniform_real_distribution<> disX(-WORLD_WIDTH/16, WORLD_WIDTH/16);
        std::uniform_real_distribution<> disY(-WORLD_HEIGTH/16, WORLD_HEIGTH/16);
        std::uniform_real_distribution<> disV(-2, 2);
        std::uniform_real_distribution<> disA(0, 2.0*M_PI);

        // Initialize particles
        particles.clear();
        particles.reserve(K_PARTICLES); //BAD choice of container if dynamic
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

            //build pairs of particles
            //if (i % 2 == 1) {
            //    particles[i].linked.push_back(&particles[i-1]);
            //    particles[i-1].linked.push_back(&particles[i]);
            //}
        }
    }

    /*
    void calculateDV_fattyAcid1(const sf::Vector2f& pos1,
                                const sf::Vector2f& pos2,
                                float orientation1,
                                float orientation2,
                                const sf::Vector2f& r,
                                float rNorm,
                                sf::Vector2f& dv,
                                float& dva)
    {
        // The interaction strength
        float strengthF = 0.2f;
        float strengthT = 0.2f;
        float teta = 0.15f;

        //Compute poles locations of the other particle
        sf::Vector2f aRVect2 = unitVectorFromAngle(orientation2+teta/2.0+M_PI/2.0);
        sf::Vector2f aRPole2 = pos2 + 5.0f * aRVect2;
        sf::Vector2f aLVect2 = unitVectorFromAngle(orientation2-teta/2.0-M_PI/2.0);
        sf::Vector2f aLPole2 = pos2 + 5.0f * aLVect2;
        float DR = norm(aRPole2 - pos1);
        float DL = norm(aLPole2 - pos1);
        float repFactor = rNorm < 5.0 ? 30.0f : 1.0f;
        if (DR < DL)
        {
            dv = repFactor*strengthF/(DR+1.0f)*(aRPole2 - pos1)/DR;
        }
        else
        {
            dv = repFactor*strengthF/(DL+1.0f)*(aLPole2 - pos1)/DL;
        }

        float deltaOrientation = orientation2 - orientation1;
        if (DR < DL)
        {
            deltaOrientation += teta;
        }
        else
        {
            deltaOrientation -= teta;
        }
        // Make sure deltaOrientation is between -pi and pi
        while (deltaOrientation > M_PI) deltaOrientation -= 2 * M_PI;
        while (deltaOrientation < -M_PI) deltaOrientation += 2 * M_PI;

        if (deltaOrientation > 0)
        {
            dva = strengthT;
        }
        else
        {
            dva = -strengthT;
        }
    }
    */

    void calculateForceAndTorque_fattyAcid1(const sf::Vector2f& pos1,
                                            const sf::Vector2f& pos2,
                                            float orientation1,
                                            float orientation2,
                                            const sf::Vector2f& r,
                                            float rNorm,
                                            sf::Vector2f& oForce,
                                            float& oTorque)
    {
        // The interaction strength
        float strengthF = 1.0f;
        float strengthT = 1.0f;
        float teta = 0.15f;

        //Compute poles locations of the other particle
        sf::Vector2f aRVect2 = unitVectorFromAngle(orientation2+teta/2.0+M_PI/2.0);
        sf::Vector2f aRPole2 = pos2 + 5.0f * aRVect2;
        sf::Vector2f aLVect2 = unitVectorFromAngle(orientation2-teta/2.0-M_PI/2.0);
        sf::Vector2f aLPole2 = pos2 + 5.0f * aLVect2;
        float DR = norm(aRPole2 - pos1);
        float DL = norm(aLPole2 - pos1);
        float repFactor = rNorm < 5.0 ? 30.0f : 1.0f;
        if (DR < DL)
        {
            float colFactor = colinearFactor(aRPole2 - pos1, unitVectorFromAngle(orientation2));
            float anisoFactor = 0.001*colFactor*colFactor*1000.0+1.0; //x times more force when "inserting"
            float distanceFactor = 1.0/((DR+1.0f)*(DR+1.0f)*(DR+1.0f));
            oForce = anisoFactor*repFactor*distanceFactor*strengthF*(aRPole2 - pos1);
        }
        else
        {
            float colFactor = colinearFactor(aLPole2 - pos1, unitVectorFromAngle(orientation2));
            float anisoFactor = 0.001*colFactor*colFactor*1000.0+1.0; //x times more force when "inserting"
            float distanceFactor = 1.0/((DL+1.0f)*(DL+1.0f)*(DL+1.0f));
            oForce = anisoFactor*repFactor*distanceFactor*strengthF*(aLPole2 - pos1);
        }

        float deltaOrientation = orientation2 - orientation1;
        if (DR < DL)
        {
            deltaOrientation += teta;
        }
        else
        {
            deltaOrientation -= teta;
        }
        // Make sure deltaOrientation is between -pi and pi
        while (deltaOrientation > M_PI) deltaOrientation -= 2 * M_PI;
        while (deltaOrientation < -M_PI) deltaOrientation += 2 * M_PI;

        if (deltaOrientation > 0)
        {
            oTorque = strengthT;
        }
        else
        {
            oTorque = -strengthT;
        }

        //Torque influence as well diminish with distance
        if (DR < DL)
        {
            oTorque = oTorque/(DR+1.0f);
        }
        else
        {
            oTorque = oTorque/(DL+1.0f);
        }

    }

    void calculateForce_linearAttraction(float forceMagnitude, const sf::Vector2f& r, float rNorm, sf::Vector2f& force) {
        force = r * (forceMagnitude / rNorm);
    }

    void calculateForce_quadraticAttraction(float forceMagnitude, const sf::Vector2f& r, float rNorm, sf::Vector2f& force) {
        force = r * (forceMagnitude / (rNorm*rNorm));
    }

    void calculateForceAndTorque_dipole(const sf::Vector2f& r,
                                        float rNorm,
                                        float orientation1,
                                        float orientation2,
                                        sf::Vector2f& oForce,
                                        float& oTorque) {

        // Construct orientation vectors
        sf::Vector2f u1 = unitVectorFromAngle(orientation1);
        sf::Vector2f u2 = unitVectorFromAngle(orientation2);

        // Compute needed quantities
        float cosTheta1 = u1.x * r.x + u1.y * r.y / rNorm;
        float cosTheta2 = u2.x * r.x + u2.y * r.y / rNorm;
        float cosTheta12 = u1.x * u2.x + u1.y * u2.y;
        //float r5 = rNorm;// * rNorm * rNorm * rNorm * rNorm;

        // Compute force
sf::Vector2f force = (3.0f * cosTheta1 * cosTheta2 - cosTheta12) * r / (rNorm);
force += (cosTheta1 * u2 + cosTheta2 * u1 - 2 * cosTheta1 * cosTheta2 * r) / (rNorm * rNorm * rNorm);

        // Compute torque
        float torque = u1.x * u2.y - u1.y * u2.x - 3 * cosTheta1 * (u1.x * r.y - u1.y * r.x) / rNorm;

        // Assign output
        oForce = force;
        oTorque = torque;
    }

/*
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
            forceMagnitude = -1000.0;

        force = r * (forceMagnitude / rNorm);

        // Calculate the torque
        // The torque magnitude depends on the relative orientation and the distance
        torque = strengthT * sinOrientation * rNorm;
    }
*/

    void step()
    {
        //const float attractionStrength = 0.05f;
        const float interactionRadius = 10.0f;
        //const float attractionMinThreshold = 5.0f;
        //const float linkingThreshold = 30.0f;

        // Calculate the force and torque on particle p due to all other particles
        for (auto& p : particles) {

            p.force = sf::Vector2f(0.0, 0.0);
            p.torque = 0.0;

            //General forces with any other particles
            for (Particle& other : particles) {
                if (&other != &p) {  // Avoid self-interaction
                    sf::Vector2f r = other.position - p.position;
                    float rNorm = norm(r);
                    if (rNorm < interactionRadius) {  // Consider only particles within the interaction radius
                        // Calculate force and torque using appropriate model
                        sf::Vector2f force(0.0, 0.0);
                        float torque = 0.0;
                        //Force model for vesicle formation
                        //calculateForceAndTorque_fattyAcid1(p.position, other.position,
                        //                                   p.orientation, other.orientation,
                        //                                   r, rNorm,
                        //                                   force, torque);

                        calculateForceAndTorque_dipole(r, rNorm, p.orientation, other.orientation, force, torque);
                        //if (p.type == 0 && other.type == 0) {
                        //    calculateForce_quadraticAttraction(40.0f, r, rNorm, force);
                        //}

                        p.force += force;
                        p.torque += torque;

                        //sf::Vector2f dv(0.0, 0.0);
                        //float dva = 0.0;
                        //calculateDV_fattyAcid1(p.position, other.position,
                        //                       p.orientation, other.orientation,
                        //                       r, rNorm,
                        //                       dv, dva);
                        //p.velocity += dv;
                        //p.angularVelocity += dva;

                        //Solid repulsion
                        //if (rNorm < 2.0*DOT_SIZE)
                        //    p.force += r * (-300.0f / rNorm);
                    }
                }
                // Floc behavior
                //p.velocity = multiply(p.velocity, 0.999);
                //p.velocity += multiply(other.velocity, 0.0001);
                //p.angularVelocity = p.angularVelocity * 0.999;
                //p.angularVelocity += other.angularVelocity * 0.0001;
            }

            //Link forces
            /*assert(p.linked.size() <= 1);
            for (Particle* otherLinked : p.linked) {
                if (otherLinked != &p) {  // Avoid self-interaction
                    sf::Vector2f r = otherLinked->position - p.position;
                    float rNorm = norm(r);
                    // Calculate force and torque using appropriate model
                    sf::Vector2f force(0.0, 0.0);
                    float torque = 0.0;
                    //Force model for vesicle formation
                    calculateForce_linearAttraction(20.0f, r, rNorm, force);
                    p.force += force;
                    p.torque += torque;
                }
            }*/

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
            //Force attract to center to incentive interactions
            p.velocity -= multiply(p.position, 0.0001);
            //Sticky dissipative space and other limits
            p.velocity = multiply(p.velocity, 0.98);
            p.angularVelocity *= 0.98;
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
    }
};

void drawModel(sf::RenderWindow& ioWindow, const Model& iModel) {
    for (const auto& p : iModel.particles)
    {
        //sf::VertexArray lines(sf::Lines, 2);
        //for (auto& other : p.linked) {
        //    lines[0].position = p.position;
        //    lines[1].position = other->position;
        //    lines[0].color = sf::Color::Green;
        //    lines[1].color = sf::Color::Green;
            //ioWindow.draw(lines);
        //}

        float tailLength = 10.0;
        sf::Vertex line[] =
        {
            sf::Vertex(p.position, sf::Color::White),
            sf::Vertex(p.position + tailLength*unitVectorFromAngle(p.orientation), sf::Color::White)
        };
        ioWindow.draw(line, 2, sf::Lines);

        //Force debug
        //sf::Vertex lineF[] =
        //{
        //    sf::Vertex(p.position, sf::Color::Red),
        //    sf::Vertex(p.position + p.force, sf::Color::Red)
        //};
        //ioWindow.draw(lineF, 2, sf::Lines);


        //Speed debug
        //sf::Vertex lineV[] =
        //{
        //    sf::Vertex(p.position, sf::Color::Green),
        //    sf::Vertex(p.position + 10.0f*p.velocity, sf::Color::Green)
        //};
        //ioWindow.draw(lineV, 2, sf::Lines);

        ioWindow.draw(p.shape);
        //ioWindow.draw(lines);

        //Pole debug
        /*
        float teta = 1.5;
        sf::Vector2f aRVect = unitVectorFromAngle(p.orientation+teta/2.0+M_PI/2.0);
        sf::Vector2f aRPole = p.position + 5.0f * aRVect;
        sf::Vector2f aLVect = unitVectorFromAngle(p.orientation-teta/2.0-M_PI/2.0);
        sf::Vector2f aLPole = p.position + 5.0f * aLVect;
        sf::CircleShape aRPoleDot(1);
        aRPoleDot.setOrigin(1, 1);
        aRPoleDot.setPosition(aRPole);
        sf::CircleShape aLPoleDot(1);
        aLPoleDot.setOrigin(1, 1);
        aLPoleDot.setPosition(aLPole);
        aRPoleDot.setFillColor(sf::Color::Red);
        aLPoleDot.setFillColor(sf::Color::Red);
        ioWindow.draw(aRPoleDot);
        ioWindow.draw(aLPoleDot);
        */
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
