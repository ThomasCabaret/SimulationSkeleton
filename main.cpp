#include <SFML/Graphics.hpp>
#include <vector>
#include <thread>
#include <random>

// Simple struct to hold particle data
struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    int type;
    std::vector<Particle*> linked;
    //Having view relating object in the model object is not ideal
    //but in SFML not rebuilding the graphical object is significantly faster
    sf::CircleShape shape;
};

static int WORLD_WIDTH = 1920;
static int WORLD_HEIGTH = 1080;
static int DOT_SIZE = 2;
static sf::Vector2f DOT_OFSET = sf::Vector2f(DOT_SIZE, DOT_SIZE);

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

    Model()
    {
        init();
    }

    void init()
    {
        // Initialize random number generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> disType(0, 15);
        std::uniform_real_distribution<> disX(-WORLD_WIDTH/2, WORLD_WIDTH/2);
        std::uniform_real_distribution<> disY(-WORLD_HEIGTH/2, WORLD_HEIGTH/2);
        std::uniform_real_distribution<> disV(-2, 2);

        // Initialize particles
        particles.clear();
        for (int i = 0; i < 1000; ++i) {
            Particle p;
            p.position = sf::Vector2f(disX(gen), disY(gen));
            p.velocity = sf::Vector2f(disV(gen), disV(gen));
            p.type = disType(gen);
            p.shape = sf::CircleShape(DOT_SIZE);
            p.shape.setPosition(p.position);
            particles.push_back(p);
        }
    }

    void step()
    {
        const float attractionStrength = 0.05f;
        const float linkingThreshold = 30.0f;

        for (auto& p : particles) {
            // Update particle position
            p.position += p.velocity;

            // Look for bound creation
            for (auto& other : particles) {
                if (&p == &other) continue;

                sf::Vector2f delta = other.position - p.position;
                float distance = std::sqrt(delta.x*delta.x + delta.y*delta.y);

                // If within linking threshold, link the particles
                if (distance < linkingThreshold && p.linked.size() < 2 && other.linked.size() < 2) {
                    p.linked.push_back(&other);
                    other.linked.push_back(&p);
                }
            }

            // Apply force
            for (auto& other : p.linked) {
                sf::Vector2f delta = other->position - p.position;
                float distance = std::sqrt(delta.x*delta.x + delta.y*delta.y);
                delta /= distance; // Normalize delta
                p.velocity += delta * attractionStrength;
            }

            // TODO Torus world is problematic to draw lines between points
            // Torus world
            //if (p.position.x > WORLD_WIDTH/2) p.position.x -= WORLD_WIDTH;
            //if (p.position.x < -WORLD_WIDTH/2) p.position.x += WORLD_WIDTH;
            //if (p.position.y > WORLD_HEIGTH/2) p.position.y -= WORLD_HEIGTH;
            //if (p.position.y < -WORLD_HEIGTH/2) p.position.y += WORLD_HEIGTH;

            // Update the particle's shape position
            p.shape.setPosition(p.position);
            p.shape.setFillColor(getColor(p.type));
        }
    }
};

void drawModel(sf::RenderWindow& ioWindow, const Model& iModel) {
    for (const auto& p : iModel.particles)
    {
        for (auto& other : p.linked) {
            sf::VertexArray lines(sf::Lines, 2);
            lines[0].position = p.position + DOT_OFSET;
            lines[1].position = other->position + DOT_OFSET;
            lines[0].color = sf::Color::White;
            lines[1].color = sf::Color::White;
            ioWindow.draw(lines);
        }

        ioWindow.draw(p.shape);
    }

}

// Main function
int main()
{
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
