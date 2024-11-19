import React, { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';
import { useDataStore } from '../store/useDataStore';
import { EnergyLevel, Transition } from '../types';

export function EnergyLevelsVisualizer() {
  const svgRef = useRef<SVGSVGElement>(null);
  const { levels, transitions, setSelectedLevel } = useDataStore();
  const [zoomDomain, setZoomDomain] = useState<{ x: [number, number]; y: [number, number] } | null>(null);
  const [deltaM, setDeltaM] = useState<'all' | 2 | 0 | -2 | 'none'>('all');
  const [fieldStrength, setFieldStrength] = useState<number>(0);
  const [showShifted, setShowShifted] = useState(false);
  const [brushing, setBrushing] = useState(false);
  const [mousePosition, setMousePosition] = useState<{ x: number, y: number } | null>(null);

  const calculateShiftedEnergy = (level: EnergyLevel, F: number): number => {
    const delta = level.electricShift;
    const C = level.magneticShift;
    return level.energy + Math.sign(C) * (-delta / 2 + Math.sqrt(Math.pow(delta / 2, 2) + Math.pow(C * F, 2)));
  };

  useEffect(() => {
    if (!svgRef.current || !levels.length) return;

    const svg = d3.select(svgRef.current);
    const svgElement = svgRef.current;
    const boundingRect = svgElement.getBoundingClientRect();
    
    const width = boundingRect.width;
    const height = boundingRect.height;
    const margin = { top: 20, right: 60, bottom: 40, left: 60 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;

    svg.attr("viewBox", `0 0 ${width} ${height}`);
    svg.selectAll("*").remove();

    // Create a clip path
    svg.append("defs")
      .append("clipPath")
      .attr("id", "plot-area")
      .append("rect")
      .attr("x", margin.left)
      .attr("y", margin.top)
      .attr("width", plotWidth)
      .attr("height", plotHeight);

    // Calculate energy values and add padding to y-axis
    const displayedLevels = levels.map(level => ({
      ...level,
      displayEnergy: showShifted ? calculateShiftedEnergy(level, fieldStrength) : level.energy
    }));

    const yExtent = d3.extent(displayedLevels, d => d.displayEnergy) as [number, number];
    const yPadding = (yExtent[1] - yExtent[0]) * 0.1; // 10% padding

    // Set scales based on zoom or full range with padding
    const xScale = d3.scaleLinear()
      .domain(zoomDomain?.x || d3.extent(displayedLevels, d => d.magneticNumber) as [number, number])
      .range([margin.left, width - margin.right]);

    const yScale = d3.scaleLinear()
      .domain(zoomDomain?.y || [yExtent[0] - yPadding, yExtent[1] + yPadding])
      .range([height - margin.bottom, margin.top]);

    // Create main plot group
    const plotGroup = svg.append("g")
      .attr("clip-path", "url(#plot-area)");

    // Add mouse position tracking
    const mousePositionGroup = svg.append("g")
      .attr("class", "mouse-position")
      .style("display", "none");

    mousePositionGroup.append("text")
      .attr("class", "position-text")
      .attr("fill", "black")
      .attr("font-size", "12px");

    // Add brush with correct coordinate system
    const brush = d3.brush()
      .extent([[margin.left, margin.top], [width - margin.right, height - margin.bottom]])
      .on("start", () => {
        setBrushing(true);
        mousePositionGroup.style("display", "none");
      })
      .on("brush", (event) => {
        if (!event.selection) return;
        const [[x0, y0], [x1, y1]] = event.selection as [[number, number], [number, number]];
        const xValue = xScale.invert(x1).toFixed(2);
        const yValue = yScale.invert(y1).toFixed(2);
        
        mousePositionGroup
          .style("display", null)
          .attr("transform", `translate(${x1 + 5}, ${y1 - 5})`);
        
        mousePositionGroup.select("text")
          .text(`2M: ${xValue}, E: ${yValue} cm⁻¹`);
      })
      .on("end", (event) => {
        setBrushing(false);
        if (!event.selection) return;

        const [[x0, y0], [x1, y1]] = event.selection as [[number, number], [number, number]];
        const newDomain = {
          x: [xScale.invert(x0), xScale.invert(x1)] as [number, number],
          y: [yScale.invert(y1), yScale.invert(y0)] as [number, number]
        };
        setZoomDomain(newDomain);
        svg.select(".brush").call(brush.move, null);
        mousePositionGroup.style("display", "none");
      });

    // Add brush group
    svg.append("g")
      .attr("class", "brush")
      .call(brush);

    // Draw axes
    svg.append("g")
      .attr("transform", `translate(0,${height - margin.bottom})`)
      .call(d3.axisBottom(xScale).ticks(5))
      .append("text")
      .attr("x", width / 2)
      .attr("y", 35)
      .attr("fill", "currentColor")
      .text("Magnetic Quantum Number (2M)");

    svg.append("g")
      .attr("transform", `translate(${margin.left},0)`)
      .call(d3.axisLeft(yScale).ticks(10))
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("x", -height / 2)
      .attr("y", -45)
      .attr("fill", "currentColor")
      .text(`Energy (cm⁻¹)${showShifted ? ' with Field' : ''}`);

    // Filter and draw transitions with correct ΔM filtering
    const filteredTransitions = transitions.filter(t => {
      if (deltaM === 'none') return false;
      if (deltaM === 'all') return true;
      
      const delta2M = t.upperState.magneticNumber - t.lowerState.magneticNumber;
      return delta2M === deltaM;
    });

    const transitionGroup = plotGroup.append("g");
    filteredTransitions.forEach(transition => {
      const upperLevel = displayedLevels.find(l => 
        l.manifold === transition.upperState.manifold &&
        l.magneticNumber === transition.upperState.magneticNumber &&
        l.stateNumber === transition.upperState.stateNumber
      );
      const lowerLevel = displayedLevels.find(l =>
        l.manifold === transition.lowerState.manifold &&
        l.magneticNumber === transition.lowerState.magneticNumber &&
        l.stateNumber === transition.lowerState.stateNumber
      );

      if (upperLevel && lowerLevel) {
        transitionGroup.append("line")
          .attr("x1", xScale(upperLevel.magneticNumber))
          .attr("y1", yScale(upperLevel.displayEnergy))
          .attr("x2", xScale(lowerLevel.magneticNumber))
          .attr("y2", yScale(lowerLevel.displayEnergy))
          .attr("stroke", "rgba(75, 85, 99, 0.2)")
          .attr("stroke-width", Math.log(transition.dipoleStrength + 1));
      }
    });

    // Draw energy levels
    const levelGroup = plotGroup.append("g");
    displayedLevels.forEach(level => {
      levelGroup.append("circle")
        .attr("cx", xScale(level.magneticNumber))
        .attr("cy", yScale(level.displayEnergy))
        .attr("r", 4)
        .attr("fill", d3.interpolateViridis(level.manifold / Math.max(...levels.map(l => l.manifold))))
        .on("mouseover", (event) => {
          if (brushing) return;
          d3.select(event.currentTarget)
            .attr("r", 6)
            .attr("stroke", "#000")
            .attr("stroke-width", 2);
          setSelectedLevel(level);
        })
        .on("mouseout", (event) => {
          if (brushing) return;
          d3.select(event.currentTarget)
            .attr("r", 4)
            .attr("stroke", "none");
          setSelectedLevel(null);
        });
    });

  }, [levels, transitions, deltaM, fieldStrength, showShifted, zoomDomain, brushing]);

  return (
    <div className="space-y-4">
      <div className="flex gap-4 items-center">
        <select
          className="px-3 py-2 border rounded-md"
          value={deltaM}
          onChange={(e) => setDeltaM(e.target.value as any)}
        >
          <option value="all">All Transitions</option>
          <option value="2">Δ(2M) = +2</option>
          <option value="0">Δ(2M) = 0</option>
          <option value="-2">Δ(2M) = -2</option>
          <option value="none">No Transitions</option>
        </select>
        
        <div className="flex items-center gap-2">
          <input
            type="checkbox"
            id="showShifted"
            checked={showShifted}
            onChange={(e) => setShowShifted(e.target.checked)}
          />
          <label htmlFor="showShifted">Show Field Shifts</label>
        </div>

        {showShifted && (
          <div className="flex items-center gap-2">
            <label htmlFor="fieldStrength">Field Strength:</label>
            <input
              type="number"
              id="fieldStrength"
              value={fieldStrength}
              onChange={(e) => setFieldStrength(Number(e.target.value))}
              className="px-3 py-2 border rounded-md w-24"
            />
          </div>
        )}

        <button
          onClick={() => setZoomDomain(null)}
          className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
        >
          Reset Zoom
        </button>
      </div>

      <div className="w-full h-[600px] bg-white rounded-lg shadow-md p-4">
        <svg
          ref={svgRef}
          className="w-full h-full"
          preserveAspectRatio="xMidYMid meet"
        />
      </div>
    </div>
  );
}