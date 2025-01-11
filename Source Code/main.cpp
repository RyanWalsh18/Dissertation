#include "body.h"
#include "quadTree.h"
#include "simulation.h"
#include <iostream>
#include <vector>
#include <SDL2/SDL.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <chrono>

using namespace nBodySimulation;

int mode = -1;
int numOfBodiesToCreate = 2;
bool defaultBodies = true;

void simMode(char *simMode)
{
    if (0 == strcmp(simMode, "brute-force"))
    {
        mode = 0;
    }

    if (0 == strcmp(simMode, "barnes-hut"))
    {
        mode = 1;
    }

    if (0 == strcmp(simMode, "fast-multipole"))
    {
        mode = 2;
    }

    if (0 == strcmp(simMode, "parallel-brute-force"))
    {
        mode = 3;
    }

    if (0 == strcmp(simMode, "parallel-barnes-hut"))
    {
        mode = 4;
    }

    if (0 == strcmp(simMode, "parallel-fast-multipole"))
    {
        mode = 5;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " < brute-force | barnes-hut | fast-multipole | parallel-brute-force | parallel-barnes-hut | parallel-fast-multipole >"
                  << " [number-of-bodies-simulated]" << std::endl;
        return 1;
    }

    if (argc == 2)
    {
        simMode(argv[1]);
        if (mode == -1)
        {
            std::cout << "Usage: " << argv[0] << " < brute-force | barnes-hut | fast-multipole | parallel-brute-force | parallel-barnes-hut | parallel-fast-multipole >"
                      << " [number-of-bodies-simulated]" << std::endl;
            return 1;
        }
    }

    if (argc == 3)
    {
        simMode(argv[1]);
        if (mode == -1)
        {
            std::cout << "Usage: " << argv[0] << " < brute-force | barnes-hut | fast-multipole | parallel-brute-force | parallel-barnes-hut | parallel-fast-multipole >"
                      << " [number-of-bodies-simulated]" << std::endl;
            return 1;
        }

        int bodiesToCreate = atoi(argv[2]);
        if (bodiesToCreate > 0)
        {
            numOfBodiesToCreate = bodiesToCreate;
            defaultBodies = false;
        }
        else
        {
            std::cout << "Usage: " << argv[0] << " < brute-force | barnes-hut | fast-multipole | parallel-brute-force | parallel-barnes-hut | parallel-fast-multipole >"
                      << " [number-of-bodies-simulated]" << std::endl;
            return 1;
        }
    }

    bool running = true;
    int WIDTH = 640;
    int HEIGHT = 640;

    int X_MID = WIDTH / 4;
    int Y_MID = HEIGHT / 4;

    int X_BOUND = WIDTH / 2 - 2;
    int Y_BOUND = HEIGHT / 2 - 2;

    SDL_Window *window = nullptr;
    SDL_Renderer *renderer = nullptr;
    SDL_Event event;
    SDL_Init(SDL_INIT_EVERYTHING);

    int FPS = 30;
    int frameDelay = 2000 / FPS;

    Uint32 frameStart;
    int frameTime;

    SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, 0, &window, &renderer);
    SDL_RenderSetScale(renderer, 2, 2);

    int finished = 0;
    // body (double bodyMass, double xPos, double yPos, double xVel, double yVel, double xAcc, double yAcc)

    // Simulation(double currentTime, double endTime, double tStep, double softeningLen, double gravConst);

    // Create bodies
    std::vector<nBodySimulation::Body> bodies;

    if (defaultBodies == false)
    {
        int numRows = std::round(std::sqrt(numOfBodiesToCreate));
        int numCols = (numOfBodiesToCreate + numRows - 1) / numRows; // Ensures all bodies are placed even if not a perfect square

        double spacingX = X_BOUND / (numCols + 1);
        double spacingY = Y_BOUND / (numRows + 1);

        for (int i = 1; i <= numRows; ++i) {
            for (int j = 1; j <= numCols; ++j) {
                if (bodies.size() >= numOfBodiesToCreate) {
                    break; // Avoid creating more bodies than needed
                }
                else{
                    double x = j * spacingX;
                    double y = i * spacingY;

                    Body body = Body(250, x, y, 0, 0); 
                    bodies.push_back(body);
                }
            }
        }
    }
    else
    {
        Body body1 = Body(10000, X_MID, Y_MID, 0, 0);
        Body body2 = Body(1, X_MID + 100, Y_MID, 0, 10.0);

        bodies.push_back(body1);
        bodies.push_back(body2);

    }

    // Set up the simulation
    nBodySimulation::Simulation simulation(0.0, 30, 0.1, 1e-6, -(1));
    // Add bodies to the simulation
    simulation.setNumOfBodies(bodies.size());
    
    auto start = std::chrono::high_resolution_clock::now();

    // Run the simulation
    while (running)
    {

        frameStart = SDL_GetTicks();

        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                running = false;
            }
        }

        // draw screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 128);

        // clear the screen
        SDL_RenderClear(renderer);
        // do drawing
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);

        for (int i = 0; i < bodies.size(); i++)
        {
            SDL_RenderDrawPoint(renderer, bodies.at(i).getX(), bodies.at(i).getY());
        }

        // show what was drawn
        SDL_RenderPresent(renderer);

        if (mode == 0)
        {
            if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
            {
                simulation.runBruteForce(bodies);
            }
            else{
                finished ++;
            }
        }

        else if (mode == 1)
        {
            if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
            {
                simulation.runBarnesHut(bodies, X_BOUND, Y_BOUND);
            }
            else{
                finished ++;
            }
        }

        else if (mode == 2)
        {
            QuadTree::precomputeBinomialCoefficients();
            QuadTree::precomputeTranslationMatrix();
            QuadTree::precomputeZPowers(24);

            if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
            {
                simulation.runFastMultipole(bodies, X_BOUND, Y_BOUND);
            }
            else{
                finished ++;
            }
        }

        else if (mode == 3)
        {
            if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
            {
                simulation.runParallelBruteForce(bodies);
            }
            else{
                finished ++;
            }
        }

        else if (mode == 4)
        {
            if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
            {
                simulation.runParallelBarnesHut(bodies, X_BOUND, Y_BOUND);
            }
            else{
                finished ++;
            }
        }

        else if (mode == 5)
        {
            QuadTree::precomputeBinomialCoefficients();
            QuadTree::precomputeTranslationMatrix();
            QuadTree::precomputeZPowers(24);

            if (simulation.getCurrentSimTime() < simulation.getSimEndTime())
            {
                simulation.runParallelFastMultipole(bodies, X_BOUND, Y_BOUND);
            }
            else{
                finished ++;
            }
        }

        if(finished == 1){
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
        }

        frameTime = SDL_GetTicks() - frameStart;
        if (frameDelay > frameTime)
        {
            SDL_Delay(frameDelay - frameTime);
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}