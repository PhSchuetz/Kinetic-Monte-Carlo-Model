# Kinetic-Monte-Carlo-Model
This is a simple Kinetic-Monte-Carlo model that describes the growth dynamics of ternary perovskite oxides (ABO<sub>3</sub>) along the (100) crystal direction during Pulsed Laser Deposition and simulates the intensity of the specular Reflection High-Energy Electron Diffraction (RHEED) signal in the step density model. 

In particular, this model is designed to describe the growth of materials like SrIrO<sub>3</sub> or SrRuO<sub>3</sub> with an unstable B-site surface termination. When grown on an A-site-terminated substrate (e.g. TiO<sub>3</sub>-terminated SrTiO<sub>3</sub>), conventionally, a A-side terminated film structure would be expected. However, if the B-cation is prone to the formation of a volatile compound upon over-oxidation (e.g. IrO<sub>3</sub> for SrIrO<sub>3</sub> or Ru<sub>2</sub>O<sub>4</sub>/RuO<sub>4</sub> for SrRu<sub>3</sub>) a termination conversion to an A-site termination can be observed during growth. Depending on the substrate temperature and oxygen partial pressure this conversion may stretch over the time needed for the deposition of one or several atomic monolayers. This model was designed to simulate the effect of this conversion onto the observed intensity oscillations of the specular RHEED signal during growth.

The basic idea of the Kinetic-Monte-Carlo simulation was adopted from the following publications:

* [P.A. Maksym, Fast Monte Carlo simulation of MBE growth](https://iopscience.iop.org/article/10.1088/0268-1242/3/6/014)

* [S. Clarke, M.R. Wilby & D.D. Vvedensky, Theory of homoepitaxy on Si(001): I. Kinetics during growth](https://www.sciencedirect.com/science/article/pii/003960289190013I#!)

* [S. Clarke & D.D. Vvedensky, Influence of surface step density on reflection high-energy-electron diffraction specular intensity during epitaxial growth](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.36.9312)

* [V.S. Achutharaman, N. Chandrasekhar, O.T. Valls, & A.M. Goldman, Origin of RHEED intensity oscillations during the growth of (Y,Dy)Ba2Cu3O7−x thin films](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.50.8122)




