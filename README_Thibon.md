Conseils d'éxecution
=
## Tests
- ***t_quad*** : Test de la quadrature de second degré (somme des poids)
- ***t_edge*** : Test de la transformée pour le point (0.2, 0) du bord 4
- ***t_triangle*** : Test de la transformée pour le point (0.2, 0.4) du triangle 4
- ***t_jac_edge*** : Test du jacobien et son déterminant pour le point (0.2, 0.4) du bord 4
- ***t_jac_triangle*** : Test du jacobien et son déterminant pour le point (0.2, 0.4) du triangle 4
- ***t_assemble_matrix*** : Test de la constitution de la matrice locale Ke pour le triangle 4 et une quadrature d'ordre 2
- ***t_assemble_vector*** : Test de la constitution du vecteur local Fe pour le triangle 4 et le bord 4 et une quadrature d'ordre 2
- ***t_Ke*** : Test du passage des matrices Ke à la matice globale K poour le maillage *square.mesh*
- ***t_F*** : Test du passage des vecteurs Fe au vecteur globale F poour le maillage *square.mesh*

## Simulations
- ***simu_pure_dirichlet*** : Simulation du problème de dirichlet pur sur le maillage *square_fine.mesh*
- ***simu_source_dirichlet*** : Simulation du problème de dirichlet avec un terme source sur le maillage *square_fine.mesh*
- ***simu_sinus_bump_dirichlet*** : Simulation du problème sinus bump sur le maillage *square.mesh* (pour avoir la **différence avec la résolution analytique**, changer *diff* en *true* dans <u>src/simu.h</u>) 
- ***simu_nuemann*** : Simulation du problème de dirichlet et de neumann sur le maillage *square_fine.mesh*
- ***simu_mug*** : Simulation du problème du mug sur le maillage *mug_1.mesh*
