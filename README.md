# asteroid_family
Auhtor: Walter S. Martins Filho
email: walter at on.br

This package is a result from my master research in Observatório Nacional (www.on.br).

Abstract:

Asteroids families are clusters of asteroids in the space of proper elements which are believed to have originated in the same parent body (Hirayama, 1918; Zappal`a, 1996). The currently method use to definy the clusters are the Hierarchical Clustering Method (HCM). This method based on the idea that asteroids that formed from the same parent body are close in the space of proper elements (Zappala & Cellino, 1990, 1994). With this, we define a dynamic family of asteroids. However, this method of identification is based on dynamic scattering of fragments after the collisional disruption of the parent body, and doesn’t infer if the asteroids really came from the same parent body. For a redefinition of the dynamics families, are used observational properties, such as color (Ivezic et al., 2002; Paolicchi & Micheli, 2008) and taxonomy (Moth ́e-Diniz et al., 2005), or the determinagion of the family age (Vokrouhlický et al., 2006). However, these methods are based on the hypothesis of the compositional homogeneity of the members of the family, and doesn’t apply to a differentiated asteroid family.

Differentiated family of asteroids is a family of asteroids that came from a differentiated parental body. The existence of metallic meteorites (McCoy et al., 2006; Elkins-Tanton et al., 2011), taxonomic diversity in asteroid families (Mothé-Diniz et al.,2005; Roig et al., 2008), and asteroids spectra compatible with differentiated achondrites (Moth ́e-Diniz & Carvano, 2005) are evidences for the existence of a differentiated asteroids families in the Main Belt. However, until now there was no real confirmation of the existence of different families (Bottke et al., 2006; Weiss & Elkins-Tanton, 2013). This leads to questioning whether the family identification methods are able to identify differences families. To test this hypothesis,we must to create a synthetic differentiated family of asteroids to be able to test the identification methods. This work propose to create a simple model that generates a different synthetic family.

The new model was based on the analytical model of (Petit & Farinella, 1993) and numerical results of Jutzi et al. (2010).

___

This package is usefull to simulate a synthethic asteroid family.

Functions:

homogeneus_family
differentiated_family
mass_distribution
min_Vej
mean_vej
mean_vej_distribution
velocity_field
catastrophic_energy
yarkovsky_dadt
yakovsky_change_unit
mag_absoluta


___

This package need:

Numpy

Scipy

Astropy

___

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
