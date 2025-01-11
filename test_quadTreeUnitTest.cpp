#include "catch.hpp"
#include "quadTree.h"
#include "body.h"
#include <complex>

using namespace nBodySimulation;

TEST_CASE("QuadTree initialisation", "[QuadTree]") {
    Quad boundary(0, 0, 100, 100);
    QuadTree qt(boundary, 4, true);


    REQUIRE(qt.getCapacity() == 4);
    REQUIRE(qt.isDivided() == false);
    REQUIRE(qt.getTotalMass() == Approx(0.0));
    REQUIRE(qt.getComX() == Approx(0.0));
    REQUIRE(qt.getComY() == Approx(0.0));
    REQUIRE(qt.isRoot() == true);
    REQUIRE(qt.getBodies().empty());
}

TEST_CASE("QuadTree body insertion and mass calculation", "[QuadTree]") {
    Quad boundary(0, 0, 200, 200);
    QuadTree qt(boundary, 4, true);
    Body body1(10, 50, 50, 0, 0);

    REQUIRE(qt.insert(body1) == true);

    REQUIRE(qt.getTotalMass() == Approx(10.0));
    REQUIRE(qt.getComX() == Approx(50.0));
    REQUIRE(qt.getComY() == Approx(50.0));
}

TEST_CASE("QuadTree subdivision functionality", "[QuadTree]") {
    Quad boundary(0, 0, 200, 200);
    QuadTree qt(boundary, 1, false); 
    Body body1(1, 25, 25, 0, 0);
    Body body2(1, 175, 175, 0, 0);

    qt.insert(body1);
    qt.insert(body2);  

    REQUIRE(qt.isDivided() == true);
    REQUIRE(qt.getNorthwest() != nullptr);
    REQUIRE(qt.getNortheast() != nullptr);
    REQUIRE(qt.getSouthwest() != nullptr);
    REQUIRE(qt.getSoutheast() != nullptr);
}

TEST_CASE("QuadTree force calculation correctness", "[QuadTree]") {
    Quad boundary(0, 0, 400, 400);
    QuadTree qt(boundary, 10, true);
    Body body1(10, 100, 100, 0, 0);
    Body body2(10, 300, 300, 0, 0);

    qt.insert(body1);
    qt.insert(body2);

    SECTION("Barnes Hut calculation between two bodies") {
        double forceX = 0.0, forceY = 0.0;
        qt.calculateForce(body1, forceX, forceY, 1.0, 6.67430e-11);

        REQUIRE(forceX != 0.0);
        REQUIRE(forceY != 0.0);
    }

    SECTION("Fast Multipole Method (FMM) force calculation") {
        double forceX = 0.0, forceY = 0.0;
        qt.calculateForcesFMM(body1, forceX, forceY, 1.0, 6.67430e-11);

        REQUIRE(forceX != 0.0);
        REQUIRE(forceY != 0.0);
    }
}
