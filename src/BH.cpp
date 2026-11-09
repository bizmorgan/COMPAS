#include "BH.h"

/*
 * Calculate the gravitational mass of the Black Hole in Msol
 *
 *
 * double CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass)
 *
 * @param   [IN]    p_BaryonicMass              Baryonic remnant mass in Msol
 * @return                                      Gravitational mass of remnant in Msol
 */
double BH::CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass) {

    double gravitationalMass = 0.0;

    switch (OPTIONS->NeutrinoMassLossAssumptionBH()) {                                                  // which assumption?

        case NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION:                                           // FIXED FRACTION
            gravitationalMass = p_BaryonicMass * (1.0 - OPTIONS->NeutrinoMassLossValueBH());
            break;

        case NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS:                                               // FIXED MASS
            gravitationalMass = p_BaryonicMass - OPTIONS->NeutrinoMassLossValueBH();
            break;
    
        default:                                                                                        // unknown prescription
            // the only way this can happen is if someone added a NEUTRINO_MASS_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose an assumption this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing assumption or, if the missing
            // assumption is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION);                         // throw error
    }

    return gravitationalMass;
}


/*
 * Calculate core collapse Supernova parameters
 *
 *
 * DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Tuple containing Luminosity, Radius and Temperature of Black Hole
 */
DBL_DBL_DBL BH::CalculateCoreCollapseSNParams_Static(const double p_Mass) {
    double luminosity  = BH::CalculateLuminosityOnPhase_Static();                                   // luminosity of BH
    double radius      = BH::CalculateRadiusOnPhase_Static(p_Mass);                                 // Schwarzschild radius (not correct for rotating BH)
    double temperature = BaseStar::CalculateTemperatureOnPhase_Static(luminosity, radius);          // temperature of BH

    return std::make_tuple(luminosity, radius, temperature);
}


/*
 * Calculate the kick given to a black hole based on users chosen assumptions about black hole kicks,
 * fallback and the magnitude of a kick drawn from a distribution
 *
 * Current options are:
 *
 *    FULL    : Black holes receive the same kicks as neutron stars
 *    REDUCED : Black holes receive the same momentum kick as a neutron star, but downweighted by the black hole mass
 *    ZERO    : Black holes receive zero natal kick
 *    FALLBACK: Black holes receive a kick downweighted by the amount of mass falling back onto them
 *
 *
 *  double ReweightSupernovaKickByMass_Static(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass)
 *
 * @param   [IN]    p_vK                        Kick magnitude that would otherwise be applied to a neutron star
 * @param   [IN]    p_FallbackFraction          Fraction of mass that falls back onto the proto-compact object
 * @param   [IN]    p_BlackHoleMass             Mass of remnant (in Msol)
 * @return                                      Kick magnitude
 */
 double BH::ReweightSupernovaKickByMass_Static(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass) {

    double vK;

    switch (OPTIONS->BlackHoleKicksMode()) {                                                            // which BH kicks mode?

        case BLACK_HOLE_KICKS_MODE::ZERO    : vK = 0.0; break;                                          // BH Kicks are set to zero regardless of BH mass or kick magnitude drawn.

        case BLACK_HOLE_KICKS_MODE::FULL    : vK = p_vK; break;                                         // BH receives full kick - no adjustment necessary

        case BLACK_HOLE_KICKS_MODE::REDUCED : vK = p_vK * NEUTRON_STAR_MASS / p_BlackHoleMass; break;   // kick is reduced by the ratio of the neutron star mass to the black hole mass

        case BLACK_HOLE_KICKS_MODE::FALLBACK: vK = p_vK * (1.0 - p_FallbackFraction); break;            // using the so-called 'fallback' mode for BH kicks
    
        default:                                                                                        // unknown mode
            // the only way this can happen is if someone added a BLACK_HOLE_KICKS_MODE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should be
            // flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_BH_KICK_MODE);                                            // throw error
    }

    return vK;
}

/*
 * Calculate the initial angular momentum of a black hole based on the angular momentum of its progenitor and the fallback fraction
 *
 * double CalculateInitialAngularMomentum_Static(const double p_ProgenitorAngularMomentum, const double p_FallbackFraction)
 *
 * @param   [IN]    p_ProgenitorAngularMomentum  Angular momentum of the progenitor star in Msol * AU^2 / yr
 * @param   [IN]    p_FallbackFraction           Fraction of mass that falls back onto the proto-compact object
 * @return                                      Initial angular momentum of the black hole in Msol * AU^2 / yr
 */

 double BH::CalculateInitialAngularMomentum_Static(
    const double p_ProgenitorAngularMomentum,
    const double p_FallbackFraction) { //this need to be replaced with my new effective fallback 

    const double fallbackFraction =
        std::max(0.0, std::min(1.0, p_FallbackFraction));

    return fallbackFraction * p_ProgenitorAngularMomentum; //simple for now, will need to apply that this is on the disk. 
}

void BH::InitialiseRemnantAngularMomentum() {

    const double progenitorAngularMomentum =
        m_SupernovaDetails.angularMomentumAtCOFormation;

    const double effectivefallbackFraction =
        GiantBranch::CalculateEffectiveFallbackFraction_Static(
            m_SupernovaDetails.preSNMass,
            m_SupernovaDetails.coreMass,
            m_SupernovaDetails.remnantMass); //biztodo check this 

    m_AngularMomentum =
        CalculateInitialAngularMomentum_Static(
            progenitorAngularMomentum,
            effectivefallbackFraction);

    // Prevent small numerical negative values
    m_AngularMomentum = std::max(m_AngularMomentum, 0.0);
}

void BH::UpdateAngularMomentumByAccretion(
    const double p_MassGain,
    const double p_SpecificAngularMomentum) {

    if (p_MassGain <= 0.0) return; //biztodo check. no accretion, no change in angular momentum. 

    const double angularMomentumGain =
        p_MassGain * p_SpecificAngularMomentum;

    SetAngularMomentum(
        AngularMomentum() + angularMomentumGain);
}

double BH::CalculateSpecificAngularMomentum_Static(
    const double p_EtaJ,
    const double p_OrbitalAngularVelocity,
    const double p_AccretionRadius) {

    return p_EtaJ *
           p_OrbitalAngularVelocity *
           p_AccretionRadius *
           p_AccretionRadius;
} //biztodo can this just go at the end? check against H&M 