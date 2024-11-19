import React from 'react';
import { useDataStore } from '../store/useDataStore';

export function LevelDetails() {
  const { selectedLevel } = useDataStore();

  if (!selectedLevel) return null;

  return (
    <div className="absolute bottom-4 right-4 p-4 bg-white rounded-lg shadow-lg max-w-md">
      <h3 className="text-lg font-semibold mb-2">Level Details</h3>
      <div className="grid grid-cols-2 gap-2 text-sm">
        <div>Manifold:</div>
        <div>{selectedLevel.manifold}</div>
        <div>Magnetic Number (2M):</div>
        <div>{selectedLevel.magneticNumber}</div>
        <div>State Number:</div>
        <div>{selectedLevel.stateNumber}</div>
        <div>Energy:</div>
        <div>{selectedLevel.energy.toFixed(6)} cm⁻¹</div>
        {selectedLevel.vibrational !== undefined && (
          <>
            <div>Vibrational Number:</div>
            <div>{selectedLevel.vibrational}</div>
          </>
        )}
        {selectedLevel.angularMomentum !== undefined && (
          <>
            <div>Angular Momentum (2J):</div>
            <div>{selectedLevel.angularMomentum}</div>
          </>
        )}
        <div>Electric Field Shift (∆):</div>
        <div>{selectedLevel.electricShift.toFixed(6)}</div>
        <div>Magnetic Field Shift (C):</div>
        <div>{selectedLevel.magneticShift.toFixed(6)}</div>
      </div>
    </div>
  );
}