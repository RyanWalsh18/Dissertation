#ifndef QUADTREE_H
#define QUADTREE_H

#include "body.h"
#include <vector>
#include <memory>
#include <complex>

namespace nBodySimulation
{
    struct MultipoleExpansion
    {
        double mass;
        std::vector<std::complex<double> > coefficients;

        MultipoleExpansion(int numCoefficients);
        MultipoleExpansion();
    };

    class Quad
    {
    private:
        double topLeftX, topLeftY;
        double bottomRightX, bottomRightY;

    public:
        Quad(double topLeftX, double topLeftY, double bottomRightX, double bottomRightY);
        bool contains(const Body &body);

        double getTopLeftX() const;
        double getTopLeftY() const;
        double getBottomRightX() const;
        double getBottomRightY() const;
    };

    class QuadTree
    {
    private:
        int capacity;
        bool divided = false;
        double totalMass = 0;
        double comX = 0, comY = 0;
        bool root = false;
        QuadTree *treeRoot;
        int level;
        double potential = 0.0;

        MultipoleExpansion outer_coefficients;
        MultipoleExpansion inner_coefficients;

        std::vector<Body> bodies;

        std::unique_ptr<Quad> boundary;

        QuadTree *parent;
        std::vector<QuadTree *> interaction_set;

        std::unique_ptr<QuadTree> northwest;
        std::unique_ptr<QuadTree> northeast;
        std::unique_ptr<QuadTree> southwest;
        std::unique_ptr<QuadTree> southeast;

        static std::vector<std::vector<double> > binomialCoefficients;
        static std::vector<std::vector<std::complex<double> > > translation_matrix;
        static std::vector<std::vector<std::complex<double> > > z_powers;

        void subdivide();
        void updateMassAndCOM(const Body &body);

        // FMM helper functions
        void compute_leaf_outer_expansion();
        void outer_shift(MultipoleExpansion &parent_expansion, const MultipoleExpansion &child_expansion, std::complex<double> &z_diff);

        void inner_shift(MultipoleExpansion &inner_expansion, const std::complex<double> &z_diff);
        void compute_interaction_set();
        void get_neighbors_at_level(std::vector<QuadTree *> &neighbors, int target_level);
        bool is_direct_neighbor(QuadTree *node);

        void add_expansions(MultipoleExpansion &target_expansion, const MultipoleExpansion &source_expansion);
        void outer_to_inner(MultipoleExpansion &inner_expansion, const MultipoleExpansion &outer_expansion, const std::complex<double> &z_diff);
        void get_direct_neighbors(std::vector<QuadTree *> &neighbors);

    public:
        QuadTree(Quad boundary, int capacity, int newLevel, QuadTree *treeRoot);
        QuadTree(Quad boundary, int capacity, bool root);

        bool insert(const Body &body);
        void calculateForce(const Body &body, double &forceX, double &forceY, double softeningLen, double gravConst);
        void build_outer_expansions();
        void build_inner_expansions();
        void calculateForcesFMM(Body &body, double &forceX, double &forceY, double softeningLen, double gravConst);

        static void precomputeBinomialCoefficients();
        static void precomputeTranslationMatrix();
        static void precomputeZPowers(int maxDegree);

        int getCapacity() const;
        bool isDivided() const;
        double getTotalMass() const;
        double getComX() const;
        double getComY() const;
        bool isRoot() const;
        int getLevel() const;
        double getPotential() const;

        const std::vector<Body>& getBodies() const;

        QuadTree* getNorthwest() const;
        QuadTree* getNortheast() const;
        QuadTree* getSouthwest() const;
        QuadTree* getSoutheast() const;
    };
}
#endif