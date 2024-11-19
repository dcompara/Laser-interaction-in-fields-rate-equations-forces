import React from 'react';
import { FileText, Users, Zap, Settings, BarChart } from 'lucide-react';

export function Overview() {
  return (
    <div className="space-y-6">
      <div className="prose max-w-none">
        <h2 className="text-2xl font-bold mb-4">Laser interaction in fields + forces</h2>
        <p>
          This application simulates the interaction of particles in various laser electric 
          and magnetic fields using rate equations.
        </p>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
        <FeatureCard
          icon={<Users className="w-6 h-6" />}
          title="Particle Configuration"
          description="Define particle types, numbers, temperatures, and initial conditions for position and velocity."
        />
        <FeatureCard
          icon={<Zap className="w-6 h-6" />}
          title="Laser & Field Settings"
          description="Configure laser parameters and external fields for precise control over the simulation environment."
        />
        <FeatureCard
          icon={<Settings className="w-6 h-6" />}
          title="Algorithm Selection"
          description="Choose between various Monte Carlo and N-body algorithms for simulation accuracy."
        />
        <FeatureCard
          icon={<BarChart className="w-6 h-6" />}
          title="Visualization Tools"
          description="Interactive Grotrian diagrams and real-time visualization of simulation results."
        />
        <FeatureCard
          icon={<FileText className="w-6 h-6" />}
          title="Data Management"
          description="Import/export configuration files and analyze simulation outputs with built-in tools."
        />
      </div>

      <div className="prose max-w-none">
        <h3 className="text-xl font-semibold mb-2">Key Features</h3>
        <ul className="list-disc pl-6">
          <li>Kinetic Monte Carlo algorithm for motion and event simulation</li>
          <li>Support for multiple particle types and laser configurations</li>
          <li>Comprehensive field control (magnetic, electric, dipolar forces)</li>
          <li>Interactive parameter configuration with real-time validation</li>
          <li>Advanced visualization tools for energy levels and transitions</li>
        </ul>
      </div>
    </div>
  );
}

function FeatureCard({ icon, title, description }: { icon: React.ReactNode, title: string, description: string }) {
  return (
    <div className="bg-white p-6 rounded-lg shadow-md">
      <div className="flex items-center mb-4">
        <div className="p-2 bg-blue-100 rounded-lg text-blue-600 mr-3">
          {icon}
        </div>
        <h3 className="text-lg font-semibold">{title}</h3>
      </div>
      <p className="text-gray-600">{description}</p>
    </div>
  );
}