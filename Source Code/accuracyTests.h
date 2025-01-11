#pragma once
void brute_force(double time_period, int numOfOrbits, double dt);
void barnes_hut(double time_period, int numOfOrbits, double dt);
void fast_multipole(double time_period, int numOfOrbits, double dt);

void parallel_brute_force(double time_period, int numOfOrbits, double dt);
void parallel_barnes_hut(double time_period, int numOfOrbits, double dt);
void parallel_fast_multipole(double time_period, int numOfOrbits, double dt);