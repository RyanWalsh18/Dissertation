#include "body.h"
#include <vector>
#include "simulation.h"
#include <iostream>
#include <cmath>
#include "quadTree.h"
#include <omp.h>

namespace nBodySimulation
{
    Simulation::Simulation()
    {
        this->numOfBodies = 0;
        this->currentSimTime = 0.0;
        this->simEndTime = 0.0;
        this->timeStep = 0.0;
        this->softeningLength = 0.0;
        this->gravConstant = 0.0;
    }
    Simulation::Simulation(double currentTime, double endTime, double tStep, double softeningLen, double gravConst)
    {
        numOfBodies = 0;
        currentSimTime = currentTime;
        simEndTime = endTime;
        timeStep = tStep;
        softeningLength = softeningLen;
        gravConstant = gravConst;
    }

    // Setter methods
    void Simulation::setNumOfBodies(int num) { numOfBodies = num; }
    void Simulation::setCurrentSimTime(double time) { currentSimTime = time; }
    void Simulation::setSimEndTime(double endTime) { simEndTime = endTime; }
    void Simulation::setTimeStep(double step) { timeStep = step; }
    void Simulation::setSofteningLength(double softening) { softeningLength = softening; }
    void Simulation::setGravConstant(double gravity) { gravConstant = gravity; }

    // Getter methods
    int Simulation::getNumOfBodies() const { return numOfBodies; }
    double Simulation::getCurrentSimTime() const { return currentSimTime; }
    double Simulation::getSimEndTime() const { return simEndTime; }
    double Simulation::getTimeStep() const { return timeStep; }
    double Simulation::getSofteningLength() const { return softeningLength; }
    double Simulation::getGravConstant() const { return gravConstant; }

    // Methods of simulation
    // FAST MULTIPOLE
    void Simulation::runFastMultipole(std::vector<Body> &bodies, int width, int height)
    {
        QuadTree quadTree(Quad(0, 0, width, height), 4, true);

        numOfBodies = bodies.size();

        for (int i = 0; i < numOfBodies; i++)
        {
            quadTree.insert(bodies.at(i));
        }

        quadTree.build_outer_expansions();

        quadTree.build_inner_expansions();

        for (int i = 0; i < numOfBodies; i++)
        {

            double forceX = 0;
            double forceY = 0;

            quadTree.calculateForcesFMM(bodies.at(i), forceX, forceY, softeningLength, gravConstant);

            forceX = -1 * (forceX);
            forceY = -1 * (forceY);

            double newXvel = bodies.at(i).getVX() + ((forceX / bodies.at(i).getMass()) * timeStep);
            double newYvel = bodies.at(i).getVY() + ((forceY / bodies.at(i).getMass()) * timeStep);

            bodies.at(i).setVX(newXvel);
            bodies.at(i).setVY(newYvel);

            double newXpos = bodies.at(i).getX() + (bodies.at(i).getVX() * timeStep);
            double newYpos = bodies.at(i).getY() + (bodies.at(i).getVY() * timeStep);

            bodies.at(i).setX(newXpos);
            bodies.at(i).setY(newYpos);
        }

        currentSimTime += timeStep;
    }

    // PARALLEL FAST MULTIPOLE
    void Simulation::runParallelFastMultipole(std::vector<Body> &bodies, int width, int height)
    {

        QuadTree quadTree(Quad(0, 0, width, height), 4, true);

        numOfBodies = bodies.size();

        for (int i = 0; i < numOfBodies; i++)
        {
            quadTree.insert(bodies.at(i));
        }

        quadTree.build_outer_expansions();
        quadTree.build_inner_expansions();
        
        std::vector<std::vector<double>> forces(numOfBodies, std::vector<double>(2, 0.0));


        int i = 0;
        #pragma omp parallel for
        for (i = 0; i < numOfBodies; i++)
        {
            quadTree.calculateForcesFMM(bodies.at(i), forces.at(i).at(0), forces.at(i).at(1), softeningLength, gravConstant);

            forces.at(i).at(0) = -1 * (forces.at(i).at(0));
            forces.at(i).at(1) = -1 * (forces.at(i).at(1));

            double newXvel = bodies.at(i).getVX() + ((forces.at(i).at(0) / bodies.at(i).getMass()) * timeStep);
            double newYvel = bodies.at(i).getVY() + ((forces.at(i).at(1) / bodies.at(i).getMass()) * timeStep);

            bodies.at(i).setVX(newXvel);
            bodies.at(i).setVY(newYvel);

            double newXpos = bodies.at(i).getX() + (bodies.at(i).getVX() * timeStep);
            double newYpos = bodies.at(i).getY() + (bodies.at(i).getVY() * timeStep);

            bodies.at(i).setX(newXpos);
            bodies.at(i).setY(newYpos);
        }

        currentSimTime += timeStep;
    }

    // BARNES HUT
    void Simulation::runBarnesHut(std::vector<Body> &bodies, int width, int height)
    {
        QuadTree quadTree(Quad(0, 0, width, height), 4, true);

        numOfBodies = bodies.size();

        for (int i = 0; i < numOfBodies; i++)
        {
            quadTree.insert(bodies.at(i));
        }

        for (int i = 0; i < numOfBodies; i++)
        {
            double forceX = 0;
            double forceY = 0;

            quadTree.calculateForce(bodies.at(i), forceX, forceY, softeningLength, gravConstant);

            forceX = -1 * (forceX);
            forceY = -1 * (forceY);

            double newXvel = bodies.at(i).getVX() + ((forceX / bodies.at(i).getMass()) * timeStep);
            double newYvel = bodies.at(i).getVY() + ((forceY / bodies.at(i).getMass()) * timeStep);

            bodies.at(i).setVX(newXvel);
            bodies.at(i).setVY(newYvel);

            double newXpos = bodies.at(i).getX() + (bodies.at(i).getVX() * timeStep);
            double newYpos = bodies.at(i).getY() + (bodies.at(i).getVY() * timeStep);

            bodies.at(i).setX(newXpos);
            bodies.at(i).setY(newYpos);
        }

        currentSimTime += timeStep;
    }

    // PARALLEL BARNES HUT
    void Simulation::runParallelBarnesHut(std::vector<Body> &bodies, int width, int height)
    {
        QuadTree quadTree(Quad(0, 0, width, height), 4, true);

        numOfBodies = bodies.size();

        for (int i = 0; i < numOfBodies; i++)
        {
            quadTree.insert(bodies.at(i));
        }

        std::vector<std::vector<double>> forces(numOfBodies, std::vector<double>(2, 0.0));

        int i = 0;
        #pragma omp parallel for
        for (i = 0; i < numOfBodies; i++)
        {
            quadTree.calculateForce(bodies.at(i), forces.at(i).at(0), forces.at(i).at(1), softeningLength, gravConstant);

            forces.at(i).at(0) = -1 * (forces.at(i).at(0));
            forces.at(i).at(1) = -1 * (forces.at(i).at(1));

            double newXvel = bodies.at(i).getVX() + ((forces.at(i).at(0) / bodies.at(i).getMass()) * timeStep);
            double newYvel = bodies.at(i).getVY() + ((forces.at(i).at(1) / bodies.at(i).getMass()) * timeStep);

            bodies.at(i).setVX(newXvel);
            bodies.at(i).setVY(newYvel);

            double newXpos = bodies.at(i).getX() + (bodies.at(i).getVX() * timeStep);
            double newYpos = bodies.at(i).getY() + (bodies.at(i).getVY() * timeStep);

            bodies.at(i).setX(newXpos);
            bodies.at(i).setY(newYpos);
        }

        currentSimTime += timeStep;
    }

    // BRUTE FORCE
    void Simulation::runBruteForce(std::vector<Body> &bodies)
    {
        numOfBodies = bodies.size();

        // Implement the simulation loop

        for (int i = 0; i < numOfBodies; i++)
        {
            for (int j = 0; j < numOfBodies; j++)
            {
                if (i != j)
                {

                    double dx = bodies.at(i).getX() - bodies.at(j).getX();
                    double dy = bodies.at(i).getY() - bodies.at(j).getY();

                    double distance = std::sqrt(dx * dx + dy * dy + softeningLength * softeningLength);

                    double forceMagnitude = (gravConstant * bodies.at(i).getMass() * bodies.at(j).getMass()) / (distance * distance);

                    double forceX = forceMagnitude * (dx / distance);
                    double forceY = forceMagnitude * (dy / distance);

                    double newXvel = bodies.at(i).getVX() + ((forceX / bodies.at(i).getMass()) * timeStep);
                    double newYvel = bodies.at(i).getVY() + ((forceY / bodies.at(i).getMass()) * timeStep);

                    bodies.at(i).setVX(newXvel);
                    bodies.at(i).setVY(newYvel);
                }
            }
        }

        for (int i = 0; i < numOfBodies; i++)
        {
            double newXpos = bodies.at(i).getX() + (bodies.at(i).getVX() * timeStep);
            double newYpos = bodies.at(i).getY() + (bodies.at(i).getVY() * timeStep);

            bodies.at(i).setX(newXpos);
            bodies.at(i).setY(newYpos);
        }

        currentSimTime += timeStep;
    }

    // PARALLEL BRUTE FORCE
    void Simulation::runParallelBruteForce(std::vector<Body> &bodies)
    {
        numOfBodies = bodies.size();

        // Implement the simulation loop

        int i;
        #pragma omp parallel for
        for (i = 0; i < numOfBodies; i++)
        {
            int j;
            for (j = 0; j < numOfBodies; j++)
            {

                if (i != j)
                {

                    double dx = bodies.at(i).getX() - bodies.at(j).getX();
                    double dy = bodies.at(i).getY() - bodies.at(j).getY();

                    double distance = std::sqrt(dx * dx + dy * dy + softeningLength * softeningLength);

                    double forceMagnitude = (gravConstant * bodies.at(i).getMass() * bodies.at(j).getMass()) / (distance * distance);

                    double forceX = forceMagnitude * (dx / distance);
                    double forceY = forceMagnitude * (dy / distance);

                    double newXvel = bodies.at(i).getVX() + ((forceX / bodies.at(i).getMass()) * timeStep);
                    double newYvel = bodies.at(i).getVY() + ((forceY / bodies.at(i).getMass()) * timeStep);

                    bodies.at(i).setVX(newXvel);
                    bodies.at(i).setVY(newYvel);
                }
            }
        }

        int n;
        #pragma omp parallel for
        for (n = 0; n < numOfBodies; n++)
        {
            double newXpos = bodies.at(n).getX() + (bodies.at(n).getVX() * timeStep);
            double newYpos = bodies.at(n).getY() + (bodies.at(n).getVY() * timeStep);

            bodies.at(n).setX(newXpos);
            bodies.at(n).setY(newYpos);
        }

        currentSimTime += timeStep;
    }
}
