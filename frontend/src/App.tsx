import React, { useState } from 'react';
import { Tabs, TabList, Tab, TabPanel } from 'react-tabs';
import { Atom, Settings, Zap, Play, Database, FileText, BarChart, Clock } from 'lucide-react';
import { Overview } from './components/Overview';
import { ParticleConfig } from './components/ParticleConfig';
import { LaserConfig } from './components/LaserConfig';
import { FieldConfig } from './components/FieldConfig';
import { AlgorithmConfig } from './components/AlgorithmConfig';
import { EnergyLevelsVisualizer } from './components/EnergyLevelsVisualizer';
import { GraphicsOutput } from './components/GraphicsOutput';
import { FileManager } from './components/FileManager';
import { ScanConfig } from './components/ScanConfig';

function App() {
  const [executablePath, setExecutablePath] = useState<string>('');

  const handleRunSimulation = () => {
    if (!executablePath) {
      alert('Please select a simulation executable in the File Management panel first.');
      return;
    }
    
    alert('Note: The simulation executable must be run locally on your computer. Web browsers cannot execute local programs directly for security reasons. Please use the desktop version of the application to run simulations.');
  };

  return (
    <div className="min-h-screen bg-gray-100">
      <header className="bg-white shadow-sm">
        <div className="max-w-7xl mx-auto px-4 py-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between">
            <div className="flex items-center">
              <Atom className="h-8 w-8 text-blue-600 mr-3" />
              <h1 className="text-2xl font-bold text-gray-900">Laser Cooling Simulation</h1>
            </div>
            <div className="flex items-center space-x-4">
              <button 
                onClick={handleRunSimulation}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 flex items-center"
              >
                <Play className="w-4 h-4 mr-2" />
                Run Simulation
              </button>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 py-6 sm:px-6 lg:px-8">
        <Tabs className="react-tabs">
          <TabList>
            <Tab><div className="flex items-center"><Database className="w-4 h-4 mr-2" />Overview</div></Tab>
            <Tab><div className="flex items-center"><Settings className="w-4 h-4 mr-2" />Particles</div></Tab>
            <Tab><div className="flex items-center"><Zap className="w-4 h-4 mr-2" />Lasers</div></Tab>
            <Tab><div className="flex items-center"><Settings className="w-4 h-4 mr-2" />Fields</div></Tab>
            <Tab><div className="flex items-center"><Settings className="w-4 h-4 mr-2" />Algorithm</div></Tab>
            <Tab><div className="flex items-center"><Zap className="w-4 h-4 mr-2" />Energy Levels</div></Tab>
            <Tab><div className="flex items-center"><BarChart className="w-4 h-4 mr-2" />Graphics & Output</div></Tab>
            <Tab><div className="flex items-center"><Clock className="w-4 h-4 mr-2" />Scan Variables</div></Tab>
            <Tab><div className="flex items-center"><FileText className="w-4 h-4 mr-2" />File Management</div></Tab>
          </TabList>

          <TabPanel><Overview /></TabPanel>
          <TabPanel><ParticleConfig /></TabPanel>
          <TabPanel><LaserConfig /></TabPanel>
          <TabPanel><FieldConfig /></TabPanel>
          <TabPanel><AlgorithmConfig /></TabPanel>
          <TabPanel><EnergyLevelsVisualizer /></TabPanel>
          <TabPanel><GraphicsOutput /></TabPanel>
          <TabPanel><ScanConfig /></TabPanel>
          <TabPanel><FileManager executablePath={executablePath} setExecutablePath={setExecutablePath} /></TabPanel>
        </Tabs>
      </main>
    </div>
  );
}

export default App;