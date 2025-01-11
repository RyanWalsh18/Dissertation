#include "catch.hpp"
#include "simulation.h"
#include "quadTree.h"

using namespace nBodySimulation;

TEST_CASE("Simulation default class constructors initialise correctly", "[simulation]") {

    Simulation sim;
    REQUIRE(sim.getNumOfBodies() == 0);
    REQUIRE(sim.getCurrentSimTime() == Approx(0.0));
    REQUIRE(sim.getSimEndTime() == Approx(0.0));
    REQUIRE(sim.getTimeStep() == Approx(0.0));
    REQUIRE(sim.getSofteningLength() == Approx(0.0));
    REQUIRE(sim.getGravConstant() == Approx(0.0));
}

TEST_CASE("Simulation parameterised class constructors initialise correctly", "[simulation]") {
    Simulation sim(1.0, 100.0, 0.1, 1e-4, 6.67430e-11);
    REQUIRE(sim.getCurrentSimTime() == Approx(1.0));
    REQUIRE(sim.getSimEndTime() == Approx(100.0));
    REQUIRE(sim.getTimeStep() == Approx(0.1));
    REQUIRE(sim.getSofteningLength() == Approx(1e-4));
    REQUIRE(sim.getGravConstant() == Approx(6.67430e-11));
}

TEST_CASE("Simulation getters and setters", "[simulation]") {
    Simulation sim(0.0, 10.0, 0.1, 1e-4, 6.67430e-11);

    SECTION("Constructor gets properties correctly") {
        REQUIRE(sim.getCurrentSimTime() == Approx(0.0));
        REQUIRE(sim.getSimEndTime() == Approx(10.0));
        REQUIRE(sim.getTimeStep() == Approx(0.1));
        REQUIRE(sim.getSofteningLength() == Approx(1e-4));
        REQUIRE(sim.getGravConstant() == Approx(6.67430e-11));
    }

    SECTION("Property setters update values correctly") {
        sim.setNumOfBodies(100);
        sim.setCurrentSimTime(1.0);
        sim.setSimEndTime(20.0);
        sim.setTimeStep(0.05);
        sim.setSofteningLength(1e-3);
        sim.setGravConstant(6.67384e-11);

        REQUIRE(sim.getNumOfBodies() == 100);
        REQUIRE(sim.getCurrentSimTime() == Approx(1.0));
        REQUIRE(sim.getSimEndTime() == Approx(20.0));
        REQUIRE(sim.getTimeStep() == Approx(0.05));
        REQUIRE(sim.getSofteningLength() == Approx(1e-3));
        REQUIRE(sim.getGravConstant() == Approx(6.67384e-11));
    }
}

TEST_CASE("Simulation methods run correctly", "[simulation]") {
    std::vector<Body> bodies = {
        Body(1.0, 100.0, 100.0, 0.0, 0.0),
        Body(1.0, 200.0, 200.0, 0.0, 0.0)
    };
    Simulation sim(0.0, 90.0, 0.1, 1e-6, -(1));
    sim.setNumOfBodies(bodies.size());

    SECTION("Brute force simulation executes without throwing exceptions") {
        REQUIRE_NOTHROW(sim.runBruteForce(bodies));
    }

    SECTION("Barnes Hut simulation executes without throwing exceptions") {
        REQUIRE_NOTHROW(sim.runBarnesHut(bodies, 800, 600));
    }

    SECTION("Fast Multipole Method simulation executes without throwing exceptions") {
        QuadTree::precomputeBinomialCoefficients();
        QuadTree::precomputeTranslationMatrix();
        QuadTree::precomputeZPowers(24);
        REQUIRE_NOTHROW(sim.runFastMultipole(bodies, 800, 600));
    }

    SECTION("Parallel brute force method execute without throwing exceptions") {
        REQUIRE_NOTHROW(sim.runParallelBruteForce(bodies));
    }

    SECTION("Parallel barnes-Hut method execute without throwing exceptions"){
        REQUIRE_NOTHROW(sim.runParallelBarnesHut(bodies, 800, 600));
    }
    
    SECTION("Parallel fast-multipole method execute without throwing exceptions"){
        QuadTree::precomputeBinomialCoefficients();
        QuadTree::precomputeTranslationMatrix();
        QuadTree::precomputeZPowers(24);
        REQUIRE_NOTHROW(sim.runParallelFastMultipole(bodies, 800, 600));
    }
}
