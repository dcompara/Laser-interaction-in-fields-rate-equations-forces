import React, { useMemo } from 'react';
import {
  useReactTable,
  getCoreRowModel,
  flexRender,
  createColumnHelper,
} from '@tanstack/react-table';
import { EnergyLevel, Transition } from '../types';
import { useDataStore } from '../store/useDataStore';

interface DataTableProps {
  type: 'levels' | 'transitions';
}

export function DataTable({ type }: DataTableProps) {
  const { levels, transitions, setLevels, setTransitions } = useDataStore();
  const columnHelper = createColumnHelper<EnergyLevel | Transition>();

  const levelColumns = useMemo(() => [
    columnHelper.accessor('manifold', {
      header: 'Manifold',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('magneticNumber', {
      header: '2M',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('boundLevel', {
      header: 'Bound',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('stateNumber', {
      header: '#',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('energy', {
      header: 'Energy (cm⁻¹)',
      cell: (info) => info.getValue().toFixed(6),
    }),
    columnHelper.accessor('electricShift', {
      header: '∆',
      cell: (info) => info.getValue().toFixed(6),
    }),
    columnHelper.accessor('magneticShift', {
      header: 'C',
      cell: (info) => info.getValue().toFixed(6),
    }),
  ], []);

  const transitionColumns = useMemo(() => [
    columnHelper.accessor('upperState.manifold', {
      header: 'Upper Manifold',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('upperState.magneticNumber', {
      header: 'Upper 2M',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('upperEnergy', {
      header: 'Upper Energy',
      cell: (info) => info.getValue().toFixed(6),
    }),
    columnHelper.accessor('lowerState.manifold', {
      header: 'Lower Manifold',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('lowerState.magneticNumber', {
      header: 'Lower 2M',
      cell: (info) => info.getValue(),
    }),
    columnHelper.accessor('lowerEnergy', {
      header: 'Lower Energy',
      cell: (info) => info.getValue().toFixed(6),
    }),
    columnHelper.accessor('dipoleStrength', {
      header: 'Dipole (Debye)',
      cell: (info) => info.getValue().toFixed(6),
    }),
  ], []);

  const table = useReactTable({
    data: type === 'levels' ? levels : transitions,
    columns: type === 'levels' ? levelColumns : transitionColumns,
    getCoreRowModel: getCoreRowModel(),
  });

  return (
    <div className="overflow-x-auto bg-white rounded-lg shadow">
      <table className="min-w-full divide-y divide-gray-200">
        <thead className="bg-gray-50">
          {table.getHeaderGroups().map(headerGroup => (
            <tr key={headerGroup.id}>
              {headerGroup.headers.map(header => (
                <th
                  key={header.id}
                  className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider"
                >
                  {flexRender(
                    header.column.columnDef.header,
                    header.getContext()
                  )}
                </th>
              ))}
            </tr>
          ))}
        </thead>
        <tbody className="bg-white divide-y divide-gray-200">
          {table.getRowModel().rows.map(row => (
            <tr key={row.id}>
              {row.getVisibleCells().map(cell => (
                <td
                  key={cell.id}
                  className="px-6 py-4 whitespace-nowrap text-sm text-gray-500"
                >
                  {flexRender(
                    cell.column.columnDef.cell,
                    cell.getContext()
                  )}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}