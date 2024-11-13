//! \file reaction_product.h
//! Data for a reaction product

#ifndef OPENMC_REACTION_PRODUCT_H
#define OPENMC_REACTION_PRODUCT_H

#include "hdf5.h"

#include "openmc/angle_energy.h"
#include "openmc/endf.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/vector.h" // for vector

#include <gsl/gsl-lite.hpp>

namespace openmc {

//==============================================================================
//! Class for secondary gamma cascade data
//==============================================================================

class Cascades{
public:
//! Secondary Gamma from HDF5
  explicit Cascades(hid_t group);

  gsl::span<const double> sample_cascade(double xi) const {
      std::size_t index = xi * lengths_.size();                //  randomly sampled index of cascade lengths
      const double* strt = energies_.data() + starts_[index];  // starting address of sampled cascade energies
      const double* end  = strt + lengths_[index];             // ending address of sampled cascade energies
      return gsl::span<const double>(strt, end);               //write for new sample 
  }                 

private:
  std::vector<std::size_t> lengths_;
  std::vector<double> energies_;
  std::vector<std::size_t> starts_;
};

//==============================================================================
//! Data for a reaction product including its yield and angle-energy
//! distributions, each of which has a given probability of occurring for a
//! given incoming energy. In general, most products only have one angle-energy
//! distribution, but for some cases (e.g., (n,2n) in certain nuclides) multiple
//! distinct distributions exist.
//==============================================================================

class ReactionProduct {
public:
  //! Emission mode for product
  enum class EmissionMode {
    prompt,  // Prompt emission of secondary particle
    delayed, // Yield represents total emission (prompt + delayed)
    total    // Delayed emission of secondary particle
  };

  using Secondary = unique_ptr<AngleEnergy>;

  //! Construct reaction product from HDF5 data
  //! \param[in] group HDF5 group containing data
  explicit ReactionProduct(hid_t group);

  //! Sample an outgoing angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(double E_in, double& E_out, double& mu, uint64_t* seed) const;

  ParticleType particle_;      //!< Particle type
  EmissionMode emission_mode_; //!< Emission mode
  double decay_rate_; //!< Decay rate (for delayed neutron precursors) in [1/s]
  unique_ptr<Function1D> yield_;      //!< Yield as a function of energy
  vector<Tabulated1D> applicability_; //!< Applicability of distribution
  vector<Secondary> distribution_;    //!< Secondary angle-energy distribution
  unique_ptr<Cascades> cascades_;     //!< Cascades from cascade generator
};

} // namespace openmc

#endif // OPENMC_REACTION_PRODUCT_H
