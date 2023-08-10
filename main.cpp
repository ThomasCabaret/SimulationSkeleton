#include "imgui.h"
#include "imgui-SFML.h"
#include <SFML/Graphics.hpp>
#include <vector>
#include <list>
#include <thread>
#include <random>
#include <cmath>
#include <iostream>
#include <cassert>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static int K_INIT_PARTICLES = 0;
static int WORLD_WIDTH = 480;//960;//1920;
static int WORLD_HEIGTH = 270;//540;//1080;
static int DOT_SIZE = 2;
static int K_NB_TYPE = 16;
static sf::Vector2f DOT_OFSET = sf::Vector2f(DOT_SIZE, DOT_SIZE);

static bool g_centerize = false;
static bool g_draw_s_interaction_radius = false;
static bool g_destroy_at_boundary = false;
static bool g_spawn_at_mouse_location= false;
static bool g_source = false;
static sf::Vector2f g_source_pos = sf::Vector2f(0, 0);
static float g_persistence = 0.5;
static float s_fps = 0;
static float g_interaction_radius = 13.0f;//10.0f;
static float g_temp_speed = 0.2f;
static float g_dt = 0.1f;
static int g_ksteps_per_frame = 10;
static float g_void_viscosity = 0.998;
static float g_void_torque_viscosity = 0.977;
static float g_containing_force = 1.0;
static float g_div_angle = 0.377;
static float g_s_f_strength = 3.21f;
static float g_s_t_strength = 0.51f;
static float g_s_f_exp_power = 3.18f;
static float g_s_f_viscosity = 1.058f;
static float g_s_t_viscosity = 1.461;
static float g_opposition_threshold = 1.426;
static float g_center_force = 0.0001f;

enum ParticleType {
    F,
    A,
    B,
    S,
    L
};

// Simple struct to hold particle data
struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    float orientation;
    //float midAngle; //debug
    float angularVelocity;
    sf::Vector2f force;
    float torque;
    ParticleType type;
    int spawnStep;
    //std::vector<Particle*> linked;
    //Having view related objects in the model object is not ideal
    //but in SFML not rebuilding the graphical object is significantly faster
    sf::CircleShape shape;
};

//////////////////////////////////////////////////////////////////////////////
float angleFromVector(sf::Vector2f v) {
    return std::atan2(v.y, v.x);
}

//////////////////////////////////////////////////////////////////////////////
sf::Vector2f unitVectorFromAngle(float iAngle) {
    return sf::Vector2f(std::cos(iAngle), std::sin(iAngle));
}

//////////////////////////////////////////////////////////////////////////////
sf::Vector2f rotateVector(const sf::Vector2f& v, float alpha) {
    float cs = std::cos(alpha);
    float sn = std::sin(alpha);

    sf::Vector2f rotatedVector;
    rotatedVector.x = v.x * cs - v.y * sn;
    rotatedVector.y = v.x * sn + v.y * cs;
    return rotatedVector;
}

float dot(const sf::Vector2f& v1, const sf::Vector2f& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

/*
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
*/

//////////////////////////////////////////////////////////////////////////////
float norm(const sf::Vector2f& vec) {
    return std::sqrt(vec.x * vec.x + vec.y * vec.y);
}

//////////////////////////////////////////////////////////////////////////////
void normalize(sf::Vector2f& vec)
{
    float normVec = norm(vec);
    if (!normVec == 0) vec /= normVec;
}

//////////////////////////////////////////////////////////////////////////////
double middleAngle(double theta1, double theta2) {
    // Make sure theta1 and theta2 are in the range of [0, 2*PI]
    theta1 = fmod(theta1, 2*M_PI);
    if (theta1 < 0)
        theta1 += 2*M_PI;

    theta2 = fmod(theta2, 2*M_PI);
    if (theta2 < 0)
        theta2 += 2*M_PI;

    // Compute difference
    double diff = theta2 - theta1;
    if (diff < -M_PI)
        diff += 2*M_PI;
    else if (diff > M_PI)
        diff -= 2*M_PI;

    // Compute middle angle
    double middle = theta1 + diff / 2.0;

    // Make sure the middle angle is in the range of [0, 2*PI]
    middle = fmod(middle, 2*M_PI);
    if (middle < 0)
        middle += 2*M_PI;

    return middle;
}

//////////////////////////////////////////////////////////////////////////////
/*sf::Color getColor(int type) {
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
}*/

sf::Color getColor(const ParticleType& iParticleType) {
    switch (iParticleType) {
    case F:
        return sf::Color::Green;
    case A:
        return sf::Color::Red;
    case B:
        return sf::Color::Magenta;
    case S:
        return sf::Color(255, 200, 0); // RGB values for yolk
    case L:
        return sf::Color(27, 227, 243); // RGB values for clear blue
    default:
        return sf::Color::White; // Return black if none of the above
    }
}


//////////////////////////////////////////////////////////////////////////////
struct Plug
{
    Grid grid;

    Plug() {
        grid.resize(PLUG_NX*PLUG_NY);
    };

    int ij2k(const sf::Vector2i& ij) const{
        return PLUG_NX*ij.y+ij.x;
    }
    sf::Vector2i k2ij(int k) const {
        int i = k % PLUG_NX;
        int j = (k-i)/PLUG_NX;
        return sf::Vector2i{ i,j };
    }

    sf::Vector2i locate(const sf::Vector2f& pos) const {
        sf::Vector2i ij;
        ij.x = floor((pos.x+.5*WORLD_WIDTH)/PLUG_DX);
        ij.y = floor((pos.y+.5*WORLD_HEIGTH)/PLUG_DY);
        if(pos.x+.5*WORLD_WIDTH < 0.)
            ij.x = 0;
        if(pos.y+.5*WORLD_HEIGTH < 0.)
            ij.y = 0;
        if(pos.x+.5*WORLD_WIDTH > WORLD_WIDTH)
            ij.x = PLUG_NX-1;
        if(pos.y+.5*WORLD_HEIGTH > WORLD_HEIGTH)
            ij.y = PLUG_NY-1;
        return ij;
    }
    Cell* getCell(const Particle& part) {
        const sf::Vector2f& pos = part.position;
        int k = ij2k(locate(pos));
        return &grid[k];
    }
    std::list<Cell*>& getNeghbourCells(const Particle& part, std::list<Cell*>& neighbour) {
        int nx = g_interaction_radius/PLUG_DX+1;
        int ny = g_interaction_radius/PLUG_DY+1;
        sf::Vector2i ij = locate(part.position);
        sf::Vector2i ijc = ij-sf::Vector2i{ nx,ny };
        for(; ijc.x <= ij.x+nx; ++ijc.x)
        {
            if(ijc.x < 0 || ijc.x >= PLUG_NX)
                continue;
            for(; ijc.y <= ij.y+ny; ++ijc.y)
            {
                if(ijc.y < 0 || ijc.y >= PLUG_NY)
                    continue;
                int k = ij2k(ijc);
                neighbour.push_back(&grid[k]);
            }
            ijc.y = ij.y-ny;
        }
        return neighbour;
    }
    void addParticle(Particle& part) {
        Cell* cell = getCell(part);
        cell->push_back(&part);
        part.cell = cell;
    }
    void removeParticle(const Particle& part) {
        Cell* cell = part.cell;
        auto it = std::find(cell->begin(),cell->end(),&part);
        cell->erase(it);
    }
    void updateCell(Particle& part) {
        const sf::Vector2f& pos = part.position;
        int k = ij2k(locate(pos));
        if(&grid[k] != part.cell)
        {
            removeParticle(part);
            addParticle(part);
        }
    }
};
//////////////////////////////////////////////////////////////////////////////

struct Model {
    std::list<Particle> particles; //BAD choice of container if dynamic
    std::random_device rd;
    std::mt19937 gen;
    int _step;

    Model() : gen(rd())
    {
        init();
    }

    void init()
    {
        _step = 0;
        // Initialize particles
        particles.clear();
        for (int i = 0; i < K_INIT_PARTICLES; ++i) {
            spawn((ParticleType)(i % 4), sf::Vector2f(0, 0), 0);
        }
    }

    void spawn(const ParticleType& iParticleType, sf::Vector2f iOrigin, int iSpawnStep)
    {
        int spawningFactor = 10;
        std::uniform_int_distribution<> disType(0, K_NB_TYPE-1);
        std::uniform_real_distribution<> disX(-WORLD_WIDTH/(2*spawningFactor), WORLD_WIDTH/(2*spawningFactor));
        std::uniform_real_distribution<> disY(-WORLD_HEIGTH/(2*spawningFactor), WORLD_HEIGTH/(2*spawningFactor));
        std::uniform_real_distribution<> disV(-0.5, 0.5);
        std::uniform_real_distribution<> disA(0, 2.0*M_PI);

        Particle p;
        p.position = sf::Vector2f(disX(gen), disY(gen));
        p.position += iOrigin;
        p.velocity = sf::Vector2f(disV(gen), disV(gen));
        p.orientation = disA(gen);
        p.angularVelocity = 0.0;
        p.type = iParticleType;
        p.shape = sf::CircleShape(DOT_SIZE);
        p.shape.setOrigin(DOT_SIZE, DOT_SIZE);
        p.shape.setPosition(p.position);
        p.spawnStep = iSpawnStep;
        particles.push_back(p);
    }

    //////////////////////////////////////////////////////////////////////////////
    void calculateForceAndTorque_polar1(Particle& p,
                                        Particle& other,
                                        const sf::Vector2f& r,
                                        float rNorm,
                                        sf::Vector2f& oForce,
                                        float& oTorque)
    {
        const sf::Vector2f& pos1 = p.position;
        const sf::Vector2f& pos2 = other.position;
        float orientation1 = p.orientation;
        float orientation2 = other.orientation;

        // The interaction strength
        float teta = orientation2 - orientation1;
        // Make sure it is between -pi and pi
        while (teta > M_PI) teta -= 2 * M_PI;
        while (teta < -M_PI) teta += 2 * M_PI;
        float phi = angleFromVector(r) - orientation1;
        float anti_phi = angleFromVector(-r) - orientation2;
        // Make sure it is between -pi and pi
        while (phi > M_PI) phi -= 2 * M_PI;
        while (phi < -M_PI) phi += 2 * M_PI;
        while (anti_phi > M_PI) anti_phi -= 2 * M_PI;
        while (anti_phi < -M_PI) anti_phi += 2 * M_PI;

        float opposition = std::max(std::abs(phi), std::abs(anti_phi));

        //Here I need a f such that:
        //f(teta, phi) = f(-teta, phi-teta+pi) //action reaction symmetry)
        //f(teta, phi) = f(-teta, -phi) //mirror symmetry
        float anisoFactor = 1.0f;//sin(teta)*sin(phi)+sin(teta)*sin(phi-teta) + 1.6f;
        float g_opposition_threshold_sp = 0.1;
        if (opposition < g_opposition_threshold-g_opposition_threshold_sp) {
            anisoFactor = -1.0f;
        } else if (opposition < g_opposition_threshold+g_opposition_threshold_sp) {
            anisoFactor = ((opposition - (g_opposition_threshold-g_opposition_threshold_sp))/(2*g_opposition_threshold_sp))*2.0-1.0; //supposed to be in -1.0 1.0
        }
        //if (((2*phi-teta+M_PI)*(2*phi-teta+M_PI)+teta*teta) < 2.0*g_div_angle*g_div_angle) anisoFactor = -1.0f;
        //2*phi-teta is an invariant angle in an interacting pair
        float rotatorFactor = sin(2*phi-teta+M_PI);//(sin(teta)*sin(phi)+sin(teta)*sin(phi-teta))/2.0f;// * M_PI;
        float distanceFactor = 1.0f / std::pow(rNorm, g_s_f_exp_power);

        oForce = rotateVector(r, rotatorFactor) * (10.0f*rotatorFactor*rotatorFactor + 1.0f) * g_s_f_strength * anisoFactor * distanceFactor;

        //Something like (phi>0)
        float midAngle = middleAngle(orientation1, orientation2);
        bool isLeft = dot(r, rotateVector(unitVectorFromAngle(midAngle), -M_PI/2.0)) > 0;

        //Debug left right
        //if(isLeft) {p.shape.setFillColor(sf::Color::Red);} else {p.shape.setFillColor(sf::Color::Blue);}

        float leftTargetOffset = (-teta-g_div_angle)/2.0;
        float righTargetOffset = (-teta+g_div_angle)/2.0;
        while (leftTargetOffset > M_PI) leftTargetOffset -= 2 * M_PI;
        while (leftTargetOffset < -M_PI) leftTargetOffset += 2 * M_PI;
        while (righTargetOffset > M_PI) righTargetOffset -= 2 * M_PI;
        while (righTargetOffset < -M_PI) righTargetOffset += 2 * M_PI;
        float dirFactor = ((isLeft && (leftTargetOffset<0)) || (!isLeft && (righTargetOffset<0))) ? 1.0f : -1.0f;
        oTorque = dirFactor*g_s_t_strength;

        //Torque influence as well diminish with distance
        oTorque = oTorque/(rNorm);

        //Mutual viscosity system (particles try to slow to the center of mass of the system)
        //TODO repulsion should be an exception!
        oForce += ((p.velocity+other.velocity)/2.0f-p.velocity)*g_s_f_viscosity;
        oTorque += ((p.angularVelocity+other.angularVelocity)/2.0f-p.angularVelocity)*g_s_t_viscosity;
    }

    //////////////////////////////////////////////////////////////////////////////
    //void calculateForce_linearAttraction(float forceMagnitude, const sf::Vector2f& r, float rNorm, sf::Vector2f& force) {
    //    force = r * (forceMagnitude / rNorm);
    //}

    //////////////////////////////////////////////////////////////////////////////
    void calculateForce_quadraticAttraction(float forceMagnitude, const sf::Vector2f& r, float rNorm, sf::Vector2f& force) {
        force = r * (forceMagnitude / (rNorm*rNorm));
    }

    //////////////////////////////////////////////////////////////////////////////
    void step()
    {
        ++_step;
        // Calculate the force and torque on particle p due to all other particles
        for (auto itp = particles.begin(); itp != particles.end();) {
            auto& p = *itp;
            p.force = sf::Vector2f(0.0, 0.0);
            p.torque = 0.0;

            //General forces with any other particles
            for (Particle& other : particles) {
                if (&other != &p) {  // Avoid self-interaction
                    sf::Vector2f r = other.position - p.position;
                    float rNorm = norm(r);
                    sf::Vector2f force(0.0, 0.0);
                    float torque = 0.0;
                    // Surfactant molecules S interaction model
                    if (p.type == ParticleType::S && other.type == ParticleType::S)
                    {
                        if (rNorm < g_interaction_radius) {  // Consider only particles within the interaction radius
                            // Calculate force and torque using appropriate model
                            //Force model for vesicle formation
                            calculateForceAndTorque_polar1(p, other,
                                                           r, rNorm,
                                                           force, torque);

                            p.force += force;
                            p.torque += torque;

                           //Solid repulsion
                           if (rNorm < 2.0*2.0*DOT_SIZE) //The first 2 is for progressive smoothing
                           {
                               float factor = std::pow(2.0*DOT_SIZE/rNorm, 9);
                                p.force += -r * factor;
                           }
                        }
                    } else if ((p.type == ParticleType::S && (other.type == ParticleType::A || other.type == ParticleType::B || other.type == ParticleType::L)) ||
                              (other.type == ParticleType::S && (p.type == ParticleType::A || p.type == ParticleType::B || p.type == ParticleType::L)))
                    {
                        // Surfactant molecules S walling model
                        if (rNorm < g_interaction_radius) {
                            // Negative for repulsion
                            //calculateForce_quadraticAttraction(-0.001, r, rNorm, force);
                            float factor = std::pow(2.0*DOT_SIZE/rNorm, 6);
                                p.force += -r * factor * (p.type != ParticleType::S ? 10.0f : 0.001f);
                            p.force += force;
                        }
                    } else if ((p.type == ParticleType::A && other.type == ParticleType::F) ||
                               (other.type == ParticleType::A && p.type == ParticleType::F))
                    {
                        //Chemical force 1: F makes B when catalysed by A
                        if (rNorm < g_interaction_radius/2.0) {
                            //Maybe TODO a commit mecahnism
                            if (p.type == ParticleType::F) {p.type = ParticleType::B; p.spawnStep = _step;}
                            if (other.type == ParticleType::F) {other.type = ParticleType::B; other.spawnStep = _step;}
                        }
                    } else if ((p.type == ParticleType::B && other.type == ParticleType::F) ||
                               (other.type == ParticleType::B && p.type == ParticleType::F))
                    {
                        //Chemical force 2: F + B makes A + S
                        if (rNorm < g_interaction_radius/2.0) {
                            p.type = ParticleType::A;
                            other.type = ParticleType::S;
                            p.spawnStep = _step;
                            other.spawnStep = _step;
                        }
                    }  else if ((p.type == ParticleType::L && other.type == ParticleType::A) ||
                               (other.type == ParticleType::A && p.type == ParticleType::L))
                    {
                        //Chemical force 3: R + A makes R + F // Test reaction limitor
                        if (rNorm < g_interaction_radius/2.0) {
                            p.type = ParticleType::F;
                            other.type = ParticleType::L;
                            p.spawnStep = _step;
                            other.spawnStep = _step;
                        }
                    }
                }
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

            if (g_destroy_at_boundary) { //p.type == ParticleType::S ||
                // Containing delete
                if ((p.position.x > WORLD_WIDTH/2) || (p.position.x < -WORLD_WIDTH/2) || (p.position.y > WORLD_HEIGTH/2) || (p.position.y < -WORLD_HEIGTH/2)) {
                    itp = particles.erase(itp);
                    continue;
                }
            } else {
                // Containing forces //could be constrained by direction of v as well
                if (p.position.x > WORLD_WIDTH/2) p.force += sf::Vector2f(-g_containing_force,0.0);
                if (p.position.x < -WORLD_WIDTH/2) p.force += sf::Vector2f(g_containing_force,0.0);
                if (p.position.y > WORLD_HEIGTH/2) p.force += sf::Vector2f(0.0,-g_containing_force);
                if (p.position.y < -WORLD_HEIGTH/2) p.force += sf::Vector2f(0.0,g_containing_force);
            }

            // Brownian motion model
            // Particle are boosted in the direction of their velocity below a given value.
            // + a rotation perturbation
            if (p.type == ParticleType::S) {
                if (norm(p.velocity) <= g_temp_speed)
                {
                    p.force += 0.01f*p.velocity / g_dt;
                }
                //std::uniform_real_distribution<> disBrownian(-0.01, 0.01);
                //p.force += sf::Vector2f(disBrownian(gen), disBrownian(gen));
            } else {
                if (norm(p.velocity) <= 1.0)
                {
                    p.force += 0.01f*p.velocity / g_dt;
                }
                std::uniform_real_distribution<> disBrownian(-0.01, 0.01);
                p.force += sf::Vector2f(disBrownian(gen), disBrownian(gen));
            }

            //Done here because it's an erase loop
            ++itp;
        }

        // Update particles from forces
        for (auto& p : particles) {
            // Calculate the acceleration and angular acceleration (assuming mass and moment of inertia = 1)
            sf::Vector2f acceleration = p.force;
            float angularAcceleration = p.torque;

            // Using naive algo
            p.velocity = p.velocity + acceleration * g_dt;
            p.angularVelocity = p.angularVelocity + angularAcceleration * g_dt;

            // Cap velocity
            float maxVelocity = 2.0;
            if (norm(p.velocity) > maxVelocity) {
                normalize(p.velocity);
                p.velocity *= maxVelocity;
            }

            // Cap angular velocity
            float maxAngularVelocity = 10.0;
            if (p.angularVelocity > maxAngularVelocity) p.angularVelocity = maxAngularVelocity;
            if (p.angularVelocity < -maxAngularVelocity) p.angularVelocity = -maxAngularVelocity;

            //Force attract to center to incentive interactions
            if (g_centerize) {
                p.velocity -= g_center_force * p.position;
            }

            //Sticky dissipative space and other limits
            p.velocity = g_void_viscosity * p.velocity;
            p.angularVelocity *= g_void_torque_viscosity;

            // Update
            p.position = p.position + p.velocity * g_dt;
            p.orientation = p.orientation + p.angularVelocity * g_dt;

            // Update the particle's shape position
            p.shape.setPosition(p.position);
            p.shape.setFillColor(getColor(p.type));
        }
    }
};

void drawModel(sf::RenderWindow& ioWindow, const sf::RectangleShape& iWorldRect, const Model& iModel) {
    for (const auto& p : iModel.particles)
    {
        // Links
        //sf::VertexArray lines(sf::Lines, 2);
        //for (auto& other : p.linked) {
        //    lines[0].position = p.position;
        //    lines[1].position = other->position;
        //    lines[0].color = sf::Color::Green;
        //    lines[1].color = sf::Color::Green;
            //ioWindow.draw(lines);
        //}

        // Tail
        if (p.type == ParticleType::S) {
            float tailLength = 10.0;
            sf::Vertex line[] =
            {
                sf::Vertex(p.position, sf::Color::White),
                sf::Vertex(p.position + tailLength*unitVectorFromAngle(p.orientation+M_PI), sf::Color::White)
            };
            ioWindow.draw(line, 2, sf::Lines);
        }

        // MidAngle debug
        //float mdLength = 30.0;
        //sf::Vertex linemd[] =
        //{
        //    sf::Vertex(p.position, sf::Color::Yellow),
        //    sf::Vertex(p.position + mdLength*unitVectorFromAngle(p.midAngle), sf::Color::Yellow)
        //};
        //ioWindow.draw(linemd, 2, sf::Lines);

        //Force debug
        //sf::Vertex lineF[] =
        //{
        //    sf::Vertex(p.position, sf::Color::Red),
        //    sf::Vertex(p.position + 10.0f*p.force, sf::Color::Red)
        //};
        //ioWindow.draw(lineF, 2, sf::Lines);

        //Speed debug
        //sf::Vertex lineV[] =
        //{
        //    sf::Vertex(p.position, sf::Color::Green),
        //    sf::Vertex(p.position + 10.0f*p.velocity, sf::Color::Green)
        //};
        //ioWindow.draw(lineV, 2, sf::Lines);


        if (iModel._step > p.spawnStep && iModel._step < p.spawnStep+g_persistence*s_fps*g_ksteps_per_frame)
        {
            float alpha = 1.0f*(iModel._step-p.spawnStep)/(g_persistence*s_fps*g_ksteps_per_frame);
            if (alpha < 0.25) {alpha = 4*alpha;}
            else if (alpha > 0.25) {alpha = 1.0-(alpha-0.25)/0.75;}
            sf::CircleShape circle(3*DOT_SIZE);  // Radius of the circle
            sf::Color color = p.shape.getFillColor();  // Get the current color
            color.a = alpha * 255;  // Set the alpha component
            circle.setFillColor(color);  // Set the color with the new alpha
            circle.setOrigin(3*DOT_SIZE, 3*DOT_SIZE);
            circle.setPosition(p.position);
            ioWindow.draw(circle);
        }
        ioWindow.draw(p.shape);
        //ioWindow.draw(lineF, 2, sf::Lines);
        //ioWindow.draw(lines);

        //Draw interaction radius
        if (p.type == ParticleType::S && g_draw_s_interaction_radius) {
            sf::CircleShape circle(g_interaction_radius);  // Radius of the circle
            circle.setFillColor(sf::Color::Transparent);  // Set the fill color to transparent
            circle.setOutlineThickness(0.1f);  // Set the outline thickness
            circle.setOutlineColor(sf::Color::Green);  // Set the outline color to green
            circle.setOrigin(g_interaction_radius, g_interaction_radius);
            circle.setPosition(p.position);
            ioWindow.draw(circle);
        }
    }
    ioWindow.draw(iWorldRect);
}

// Main function
int main()
{
    // Create the main window
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Particle system");
    const bool initOK = ImGui::SFML::Init(window);
    if(!initOK){
        std::cerr << "Problem" << std::endl;
        return 1;
    }
    ImGui::GetStyle().ScaleAllSizes(2.0f);
    ImGui::GetIO().FontGlobalScale = 2.0f;

    // Model
    Model myModel;

    // Create a view with the same size as the window
    sf::View view(sf::FloatRect(-WORLD_WIDTH/2, -WORLD_HEIGTH/2, WORLD_WIDTH, WORLD_HEIGTH));

    // Variables to store the state of mouse dragging
    bool isDragging = false;
    sf::Vector2f lastMousePos;
    sf::Vector2i mouseDownPxPos;

    // Create a sf::RectangleShape with the given world size
    sf::RectangleShape worldRect(sf::Vector2f(WORLD_WIDTH, WORLD_HEIGTH));
    worldRect.setOutlineThickness(4.0f); // Adjust as needed
    worldRect.setOutlineColor(sf::Color::Red);
    worldRect.setFillColor(sf::Color::Transparent);
    worldRect.setOrigin(WORLD_WIDTH / 2.0f, WORLD_HEIGTH / 2.0f);
    worldRect.setPosition(0.0f, 0.0f);


    // main loop
    sf::Clock deltaClock;
    while (window.isOpen()) {
        // Handle events
        sf::Event event;
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(event);
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
                if (!ImGui::GetIO().WantCaptureMouse) {
                    if (event.mouseButton.button == sf::Mouse::Left) {
                        isDragging = true;
                        mouseDownPxPos = sf::Mouse::getPosition(window);
                        lastMousePos = window.mapPixelToCoords(mouseDownPxPos);
                    }
                    break;
                }
            case sf::Event::MouseButtonReleased:
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    if (!ImGui::GetIO().WantCaptureMouse) {
                        const sf::Vector2i mouseUpPxPos = sf::Mouse::getPosition(window);
                        const sf::Vector2i delta = mouseDownPxPos - mouseUpPxPos;
                        int CLICK_THRESHOLD = 2;
                        if (abs(mouseDownPxPos.x - mouseUpPxPos.x) < CLICK_THRESHOLD && abs(mouseDownPxPos.y - mouseUpPxPos.y) < CLICK_THRESHOLD) {
                            g_source = !g_source;
                            g_source_pos = window.mapPixelToCoords(mouseUpPxPos);
                        }
                    }
                    isDragging = false;
                }
                break;
            case sf::Event::KeyPressed:
                sf::Vector2f mousePxPosForSpawn;
                if(g_spawn_at_mouse_location){
                    //get current mouse position to spawn particule where mouse is
                    mousePxPosForSpawn = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                }else{
                    mousePxPosForSpawn = sf::Vector2f(0, 0);
                }
                if (event.key.code == sf::Keyboard::R)
                {
                    for(Particle& p : myModel.particles)
                        myModel.plug.removeParticle(p);
                    myModel.particles.clear();
                }
                if (event.key.code == sf::Keyboard::C)
                    g_centerize = !g_centerize;
                if (event.key.code == sf::Keyboard::S)
                    myModel.spawn(ParticleType::S, mousePxPosForSpawn, myModel._step);
                if (event.key.code == sf::Keyboard::F)
                    myModel.spawn(ParticleType::F, mousePxPosForSpawn, myModel._step);
                if (event.key.code == sf::Keyboard::A)
                    myModel.spawn(ParticleType::A, mousePxPosForSpawn, myModel._step);
                if (event.key.code == sf::Keyboard::B)
                    myModel.spawn(ParticleType::B, mousePxPosForSpawn, myModel._step);
                if (event.key.code == sf::Keyboard::L)
                    myModel.spawn(ParticleType::L, mousePxPosForSpawn, myModel._step);
                break;
            }
        }
        ImGui::SFML::Update(window, deltaClock.restart());

        //GUI
        ImGui::Begin("Demo window");
        s_fps = 1.0f / ImGui::GetIO().DeltaTime;
        ImGui::Text("FPS: %.1f", s_fps);
        ImGui::Text("Population: %d", myModel.particles.size());
        if (ImGui::InputInt("Steps per frame", &g_ksteps_per_frame)) {
            if (g_ksteps_per_frame < 1) {
                g_ksteps_per_frame = 1;
            }
        }
        ImGui::SliderFloat("dt", &g_dt, 0.0f, 1.0f);
        ImGui::SliderFloat("Containing force", &g_containing_force, 0.0f, 2.0f);
        ImGui::SliderFloat("Centerize", &g_center_force, 0.0f, 0.001f);
        ImGui::SliderFloat("Brownian speed", &g_temp_speed, 0.0f, 1.0f);
        ImGui::SliderFloat("Void viscosity", &g_void_viscosity, 0.99f, 1.0f);
        ImGui::SliderFloat("Void torque viscosity", &g_void_torque_viscosity, 0.9f, 1.0f);
        ImGui::SliderFloat("S interaction radius", &g_interaction_radius, 4.0f, 20.0f);
        ImGui::SliderFloat("S force strength", &g_s_f_strength, 0.0f, 5.0f);
        ImGui::SliderFloat("S force exp power", &g_s_f_exp_power, 1.0f, 5.0f);
        ImGui::SliderFloat("S torque strength", &g_s_t_strength, 0.0f, 2.0f);
        ImGui::SliderFloat("S div angle", &g_div_angle, 0.0, 0.4);
        ImGui::SliderFloat("S force viscosity", &g_s_f_viscosity, 0.0, 4.0f);
        ImGui::SliderFloat("S torque viscosity", &g_s_t_viscosity, 0.0, 4.0f);
        ImGui::SliderFloat("S opposition threshold", &g_opposition_threshold, 0.0, 7.0f);
        ImGui::Checkbox("Destroy at boundary", &g_destroy_at_boundary);
        ImGui::Checkbox("Draw S Interaction Radius", &g_draw_s_interaction_radius);
        ImGui::Checkbox("Spawn particules at mouse location", &g_spawn_at_mouse_location);
        ImGui::Text("S,F,A,B key to spawn particles");
        ImGui::Text("C key to center");
        ImGui::Text("R key to reset");
        ImGui::End();

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

        // Spawning
        if (g_source)
        {
            myModel.spawn(ParticleType::F, g_source_pos, myModel._step);
        }

        // Model update
        for (int i=0; i<g_ksteps_per_frame; ++i)
            myModel.step();

        // Clear screen
        window.clear();

        // Draw particles
        drawModel(window, worldRect, myModel);
        ImGui::SFML::Render(window);

        // Update the window
        window.display();
    }

    ImGui::SFML::Shutdown();

    return 0;
}
