#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "body.h"
#include <SDL2/SDL.h>

namespace nBodySimulation
{
    class Simulation
    {
    private:
        int numOfBodies;
        double currentSimTime;
        double simEndTime;
        double timeStep;
        double softeningLength;
        double gravConstant;

    public:
        // Constructors
        Simulation();

        Simulation(double currentTime, double endTime, double tStep, double softeningLen, double gravConst);

        // Setter methods
        void setNumOfBodies(int num);
        void setCurrentSimTime(double currentTime);
        void setSimEndTime(double endTime);
        void setTimeStep(double tStep);
        void setSofteningLength(double softeningLen);
        void setGravConstant(double gravConst);

        // Getter methods
        int getNumOfBodies() const;
        double getCurrentSimTime() const;
        double getSimEndTime() const;
        double getTimeStep() const;
        double getSofteningLength() const;
        double getGravConstant() const;

        void runBruteForce(std::vector<Body> &bodies);
        void runBarnesHut(std::vector<Body> &bodies, int width, int height);
        void runFastMultipole(std::vector<Body> &bodies, int width, int height);

        void runParallelBruteForce(std::vector<Body> &bodies);
        void runParallelBarnesHut(std::vector<Body> &bodies, int width, int height);
        void runParallelFastMultipole(std::vector<Body> &bodies, int width, int height);
    };
}
#endif
