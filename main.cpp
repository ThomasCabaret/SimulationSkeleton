#include <SFML/Graphics.hpp>
#include <vector>
#include <thread>
#include <random>

// Simple struct to hold particle data
struct Particle {
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::CircleShape shape;
};

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
        std::uniform_real_distribution<> dis(0, 1920);
        std::uniform_real_distribution<> disV(-2, 2);

        // Initialize particles
        particles.clear();
        for (int i = 0; i < 1000; ++i) {
            Particle p;
            p.position = sf::Vector2f(dis(gen), dis(gen));
            p.velocity = sf::Vector2f(disV(gen), disV(gen));
            p.shape = sf::CircleShape(2);
            p.shape.setPosition(p.position);
            particles.push_back(p);
        }
    }

    void step()
    {
        for (auto& p : particles) {
            // Update particle position
            p.position += p.velocity;

            // If the particle hits the boundaries, reverse its velocity
            if (p.position.x < 0 || p.position.x > 1920) {
                p.velocity.x *= -1;
            }
            if (p.position.y < 0 || p.position.y > 1080) {
                p.velocity.y *= -1;
            }

            // Update the particle's shape position
            p.shape.setPosition(p.position);
        }
    }
};

void drawModel(sf::RenderWindow& ioWindow, const Model& iModel) {
    for (const auto& p : iModel.particles)
        ioWindow.draw(p.shape);
}

// Main function
int main()
{
    // Create the main window
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Particle system");

    // Model
    Model myModel;

    // Create a view with the same size as the window
    sf::View view(sf::FloatRect(0, 0, 1920, 1080));

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
