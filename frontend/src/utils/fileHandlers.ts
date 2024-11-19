import { EnergyLevel, Transition } from '../types';

export function parseLevelsFile(content: string): EnergyLevel[] {
  const lines = content.trim().split('\n');
  return lines.map(line => {
    const [
      manifold, magneticNumber, boundLevel, stateNumber,
      population, vibrational, angularMomentum, energy,
      electricShift, magneticShift
    ] = line.trim().split(/\s+/).map(Number);

    return {
      manifold,
      magneticNumber,
      boundLevel,
      stateNumber,
      population,
      vibrational,
      angularMomentum,
      energy,
      electricShift,
      magneticShift
    };
  });
}

export function parseTransitionsFile(content: string): Transition[] {
  const lines = content.trim().split('\n');
  return lines.map(line => {
    const [
      upperManifold, upperM, upperBound, upperState,
      lowerManifold, lowerM, lowerBound, lowerState,
      upperEnergy, lowerEnergy, dipoleStrength
    ] = line.trim().split(/\s+/).map(Number);

    return {
      upperState: {
        manifold: upperManifold,
        magneticNumber: upperM,
        boundLevel: upperBound,
        stateNumber: upperState
      },
      lowerState: {
        manifold: lowerManifold,
        magneticNumber: lowerM,
        boundLevel: lowerBound,
        stateNumber: lowerState
      },
      upperEnergy,
      lowerEnergy,
      dipoleStrength
    };
  });
}