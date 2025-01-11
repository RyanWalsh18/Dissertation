#include "catch.hpp"
#include "body.h"

using namespace nBodySimulation;

TEST_CASE("Body class default constructor initialises correctly", "[body]"){
    Body body;

    REQUIRE(body.getMass() == 0.0);
    REQUIRE(body.getX() == 0.0);
    REQUIRE(body.getY() == 0.0);
    REQUIRE(body.getVX() == 0.0);
    REQUIRE(body.getVY() == 0.0);
    REQUIRE(body.getAX() == 0.0);
    REQUIRE(body.getAY() == 0.0);
}

TEST_CASE("Body class parameterized constructors initialise correctly", "[body]") {

    SECTION("Constructor with mass, position, and velocity initialises correctly") {
        Body body(1.0, 2.0, 3.0, 4.0, 5.0);
        REQUIRE(body.getMass() == 1.0);
        REQUIRE(body.getX() == 2.0);
        REQUIRE(body.getY() == 3.0);
        REQUIRE(body.getVX() == 4.0);
        REQUIRE(body.getVY() == 5.0);
        REQUIRE(body.getAX() == 0.0); 
        REQUIRE(body.getAY() == 0.0); 
    }

    SECTION("Full constructor with all parameters initialises correctly") {
        Body body(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
        REQUIRE(body.getMass() == 1.0);
        REQUIRE(body.getX() == 2.0);
        REQUIRE(body.getY() == 3.0);
        REQUIRE(body.getVX() == 4.0);
        REQUIRE(body.getVY() == 5.0);
        REQUIRE(body.getAX() == 6.0);
        REQUIRE(body.getAY() == 7.0);
    }
}

TEST_CASE("Body class setters and getters work correctly", "[body]") {
    Body body;
    body.setMass(10.0);
    body.setX(20.0);
    body.setY(30.0);
    body.setVX(40.0);
    body.setVY(50.0);
    body.setAX(60.0);
    body.setAY(70.0);

    SECTION("Setters correctly modify Body properties") {
        REQUIRE(body.getMass() == 10.0);
        REQUIRE(body.getX() == 20.0);
        REQUIRE(body.getY() == 30.0);
        REQUIRE(body.getVX() == 40.0);
        REQUIRE(body.getVY() == 50.0);
        REQUIRE(body.getAX() == 60.0);
        REQUIRE(body.getAY() == 70.0);
    }

    SECTION("setPosition correctly modifies x and y coordinates") {
        body.setPosition(100.0, 200.0);
        REQUIRE(body.getX() == 100.0);
        REQUIRE(body.getY() == 200.0);
    }

    SECTION("setVelocity correctly modifies vx and vy") {
        body.setVelocity(300.0, 400.0);
        REQUIRE(body.getVX() == 300.0);
        REQUIRE(body.getVY() == 400.0);
    }

    SECTION("setAcceleration correctly modifies ax and ay") {
        body.setAcceleration(500.0, 600.0);
        REQUIRE(body.getAX() == 500.0);
        REQUIRE(body.getAY() == 600.0);
    }
}

TEST_CASE("Comparison operators for Body class work correctly", "[body]") {
    Body body1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    Body body2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    Body body3(1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1); 

    SECTION("Bodies with the same values are equal") {
        REQUIRE(body1 == body2);
    }

    SECTION("Bodies with different values are not equal") {
        REQUIRE(body1 != body3);
    }
}

TEST_CASE("Body class handles negative and extreme values correctly", "[body]") {
    Body body;

    SECTION("Handling negative values") {
        body.setMass(-1.0);
        body.setPosition(-100.0, -200.0);
        body.setVelocity(-300.0, -400.0);
        body.setAcceleration(-500.0, -600.0);

        REQUIRE(body.getMass() == -1.0);
        REQUIRE(body.getX() == -100.0);
        REQUIRE(body.getY() == -200.0);
        REQUIRE(body.getVX() == -300.0);
        REQUIRE(body.getVY() == -400.0);
        REQUIRE(body.getAX() == -500.0);
        REQUIRE(body.getAY() == -600.0);
    }

    SECTION("Handling extreme values") {
        double largeValue = std::numeric_limits<double>::max();
        double smallValue = std::numeric_limits<double>::min();

        body.setMass(largeValue);
        body.setPosition(smallValue, largeValue);
        body.setVelocity(largeValue, smallValue);
        body.setAcceleration(smallValue, largeValue);

        REQUIRE(body.getMass() == largeValue);
        REQUIRE(body.getX() == smallValue);
        REQUIRE(body.getY() == largeValue);
        REQUIRE(body.getVX() == largeValue);
        REQUIRE(body.getVY() == smallValue);
        REQUIRE(body.getAX() == smallValue);
        REQUIRE(body.getAY() == largeValue);
    }
}
