#include <stdint.h>
void solve(uint16_t element_count, uint16_t element_nodes,
           uint16_t integration_points, uint16_t node_count,
           uint16_t plane_stress, uint16_t conditions_count,
           uint16_t reactions_count, double nu, double E, double bforce,
           const double *nodes, const uint16_t *element,
           const uint16_t *reactions, const uint16_t *displacement_nodes,
           const double *displacement_given, double *displacement_result,
           double *reaction_result);
