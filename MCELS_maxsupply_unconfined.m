function out = MCTC_isolated_50pc_sphere_unconfined(argin)
%
% MCTC_fluid_2D_90pc_sphere.m
% input arguments: 
% 1. inlet pressure
% 2. water diffusivity
% 3. MCTC diffusivity
% 4. consumption rate
% 5. Damköhler (with MCTC diffusion)
% 6. diffusivity ratio (R_d)
% 7. Péclet number (fluid channel)
%
% Model exported on Nov 6 2023, 10:05 by COMSOL 6.1.0.252.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\Willy\Work\Comsol_simulations\MCTC_reservoir_rigid\2D');

model.label('MCTC_fluid_2D_50pc_sphere_unconfined_m.mph');

model.param.set('L', '3[mm]', 'device channel length');
model.param.set('a', '0.5[mm]', 'MCTC radius');
model.param.set('b', '0 [mm]', 'height of MCTC center above bottom');
model.param.set('H', 'a+b+0.3[mm]', 'Microfluidic device height');
model.param.set('th', ['280[' native2unicode(hex2dec({'00' 'b0'}), 'unicode') ']'], 'MCTC angle');
model.param.set('k_o', '1e-16[m^2]', 'MCTC permeability');
model.param.set('eps_o', '0.9', 'MCTC porosity');
model.param.set('D_ext', '1e-9[m^2/s]', ['Molecular diffusivity in water at 37' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C']);
model.param.set('D_int', strcat(num2str(argin(1)), '[m^2/s]'), 'Diffusivity in MCTC');
model.param.set('R0', strcat(num2str(argin(2)), '[mol/m^3/s]'), 'Max. consumption rate');
model.param.set('p_0', '1e-8[Pa]', 'Inlet pressure');
model.param.set('U_0', '1[mm/s]', 'Inlet velocity');
model.param.set('D_bmp4', '9.3e-11[m^2/s]', 'Bmp4 diffusivity (water)');
model.param.set('Pr_min', '1e-11[mol/m^3/s]', 'min. production rate');
model.param.set('Pr_max', '1e-9[mol/m^3/s]', 'max. production rate');
model.param.set('Pr_C', '1e-10[mol/m^3]', '"concentration constant" of concentration-production feedback');
model.param.set('Phi_bmp4', '1e-11[mol/m^3/s]', 'Bmp4 consumption rate (assumed)');
model.param.set('MU', '0.0025[Pa*s]', 'dynamic viscosity');
model.param.set('K_m', '4.63e-5[mol/m^3]', 'Michaelis-Menten constant for O2 (concentration of half consumption)');
model.param.set('c_0', '0.2[mol/m^3]');

model.component.create('comp1', true);
model.component('comp1').label('MCTC 2D');

model.component('comp1').variable.create('var7');
model.component('comp1').variable('var7').set('r_M', 'sqrt(x^2+(y-b)^2)', 'radius in MCTC centered coordinates');
model.component('comp1').variable('var7').set('Diff', '(r_M<a)*D_int+(r_M>=a)*D_ext', 'diffusivity');
model.component('comp1').variable('var7').label('diffusivity/velocity entire model');

% geometry
model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'-L' '0'});
model.component('comp1').geom('geom1').feature('r1').set('size', {'2*L' 'H'});
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('pos', {'0' 'b'});
model.component('comp1').geom('geom1').feature('c1').set('rot', '0');
model.component('comp1').geom('geom1').feature('c1').set('type', 'curve');
model.component('comp1').geom('geom1').feature('c1').set('r', 'a');
model.component('comp1').geom('geom1').feature('c1').set('angle', '180');

model.component('comp1').geom('geom1').run;

% material
model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material('mat3').label('Culture medium');
model.component('comp1').material('mat3').propertyGroup('def').set('density', '1030');
model.component('comp1').material('mat3').propertyGroup('def').set('dynamicviscosity', 'MU');
model.component('comp1').material('mat3').materialType('nonSolid');

% laminar flow
model.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
model.component('comp1').physics('spf').identifier('spf');
model.component('comp1').physics('spf').field('velocity').field('v2');
model.component('comp1').physics('spf').field('velocity').component({'v2x' 'v2y' 'v2z'});
model.component('comp1').physics('spf').field('pressure').field('p2');
model.component('comp1').physics('spf').field('turbulentkineticenergy').field('k2');
model.component('comp1').physics('spf').field('turbulentdissipationrate').field('ep2');
model.component('comp1').physics('spf').field('specificdissipationrate').field('om2');
model.component('comp1').physics('spf').field('reciprocallength').field('G2');
model.component('comp1').physics('spf').field('correctedvelocity').field('uc2');
model.component('comp1').physics('spf').field('correctedvelocity').component({'uc2x' 'uc2y' 'uc2z'});
model.component('comp1').physics('spf').field('correctedpressure').field('pc2');
model.component('comp1').physics('spf').field('turbulentkinematicviscosity').field('nutilde2');
model.component('comp1').physics('spf').field('dimensionless1').field('yPlus2');
model.component('comp1').physics('spf').field('dimensionless2').field('uPlus2');
model.component('comp1').physics('spf').field('dimensionless3').field('zeta2');
model.component('comp1').physics('spf').field('dimensionless4').field('alpha2');
model.component('comp1').physics('spf').field('dimensionless5').field('gamma2');
model.component('comp1').physics('spf').selection.set([1]);

model.component('comp1').physics('spf').create('inl1', 'InletBoundary', 1);
model.component('comp1').physics('spf').feature('inl1').selection.set([1]);
model.component('comp1').physics('spf').feature('inl1').set('U0in', 'p_0/4/MU/L*(H^2-y^2)');

model.component('comp1').physics('spf').create('out1', 'OutletBoundary', 1);
model.component('comp1').physics('spf').feature('out1').selection.set([7]);
model.component('comp1').physics('spf').feature('out1').set('NormalFlow', true);
model.component('comp1').physics('spf').feature('out1').set('SuppressBackflow', false);

model.component('comp1').physics('spf').create('wallbc2', 'WallBC', 1);
model.component('comp1').physics('spf').feature('wallbc2').selection.set([8 9]);
model.component('comp1').physics('spf').feature('wallbc2').set('BoundaryCondition', 'LeakingWall');
model.component('comp1').physics('spf').feature('wallbc2').set('uleak', {'dl.u'; 'dl.v'; '0'}); % HERE: PROBLEM?
model.component('comp1').physics('spf').feature('wallbc2').set('UseViscousSlip', true);
model.component('comp1').physics('spf').feature('wallbc2').label('Wall interface');

model.component('comp1').physics('spf').create('sym1', 'Symmetry', 1);
model.component('comp1').physics('spf').feature('sym1').selection.set([2 6]);

model.component('comp1').physics('spf').prop('ShapeProperty').set('order_fluid', 2);
model.component('comp1').physics('spf').prop('PhysicalModelProperty').set('EnablePorousMediaDomains', true);

% Darcy's law
model.component('comp1').physics.create('dl', 'PorousMediaFlowDarcy', 'geom1');
model.component('comp1').physics('dl').identifier('dl');
model.component('comp1').physics('dl').field('pressure').field('p2d');
model.component('comp1').physics('dl').selection.set([2]);

model.component('comp1').physics('dl').feature('porous1').feature('fluid1').set('mu_mat', 'userdef');
model.component('comp1').physics('dl').feature('porous1').feature('fluid1').set('mu', 'MU');
model.component('comp1').physics('dl').feature('porous1').feature('pm1').set('epsilon_mat', 'userdef');
model.component('comp1').physics('dl').feature('porous1').feature('pm1').set('epsilon', 'eps_o');
model.component('comp1').physics('dl').feature('porous1').feature('pm1').set('kappa_mat', 'userdef');
model.component('comp1').physics('dl').feature('porous1').feature('pm1').set('kappa', {'k_o'; '0'; '0'; '0'; 'k_o'; '0'; '0'; '0'; 'k_o'});

model.component('comp1').physics('dl').create('pr1', 'Pressure', 1);
model.component('comp1').physics('dl').feature('pr1').selection.set([8 9]);
model.component('comp1').physics('dl').feature('pr1').set('p0', 'p2');

model.component('comp1').physics('dl').create('sym1', 'Symmetry', 1);
model.component('comp1').physics('dl').feature('sym1').selection.set([4 5]);

% transport of diluted species
model.component('comp1').physics.create('tds', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds').identifier('tds');
model.component('comp1').physics('tds').field('concentration').field('c2');
model.component('comp1').physics('tds').field('concentration').component({'c2'});
model.component('comp1').physics('tds').selection.set([2]);
model.component('comp1').physics('tds').prop('TransportMechanism').set('MassTransferInPorousMedia', true);

model.component('comp1').physics('tds').feature('cdm1').set('u_src', 'root.comp1.dl.u');
model.component('comp1').physics('tds').feature('cdm1').set('DiffusionCoefficientSource', 'chem');
model.component('comp1').physics('tds').feature('cdm1').set('Dchem_c2', {'D_int'; '0'; '0'; '0'; 'D_int'; '0'; '0'; '0'; 'D_int'});

model.component('comp1').physics('tds').feature('init1').set('initc', 'c_0');

% model.component('comp1').physics('tds').create('cdm2', 'ConvectionDiffusionMigration', 2);
% model.component('comp1').physics('tds').feature('cdm2').selection.set([2]);
% model.component('comp1').physics('tds').feature('cdm2').set('u_src', 'root.comp1.dl.u');
% model.component('comp1').physics('tds').feature('cdm2').set('DiffusionCoefficientSource', 'chem');
% model.component('comp1').physics('tds').feature('cdm2').set('Dchem_c2', {'D_int'; '0'; '0'; '0'; 'D_int'; '0'; '0'; '0'; 'D_int'});

model.component('comp1').physics('tds').create('conc1', 'Concentration', 1);
model.component('comp1').physics('tds').feature('conc1').selection.set([8 9]);
model.component('comp1').physics('tds').feature('conc1').set('species', true);
model.component('comp1').physics('tds').feature('conc1').set('c0', 'c_0');

model.component('comp1').physics('tds').create('reac1', 'Reactions', 2);
model.component('comp1').physics('tds').feature('reac1').selection.set([2]);
model.component('comp1').physics('tds').feature('reac1').set('R_c2', '-R0*c2/(c2+K_m)*(c2>0)');
model.component('comp1').physics('tds').feature('reac1').label('Consumption');

model.component('comp1').physics('tds').create('sym1', 'Symmetry', 1);
model.component('comp1').physics('tds').feature('sym1').selection.set([4 5]);

% mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.025);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.001);
model.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1.15);

model.component('comp1').mesh('mesh1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('size1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('size1').label('Size MCTC');
model.component('comp1').mesh('mesh1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size1').set('hmax', 0.025/sqrt(2));
model.component('comp1').mesh('mesh1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size1').set('hmin', 6.0E-4);
model.component('comp1').mesh('mesh1').feature('size1').set('hminactive', false);

model.component('comp1').mesh('mesh1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').feature('size2').selection.geom('geom1', 1);
model.component('comp1').mesh('mesh1').feature('size2').selection.set([8 9]);
model.component('comp1').mesh('mesh1').feature('size2').label('Size interface');
model.component('comp1').mesh('mesh1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size2').set('hmax', 0.00625);
model.component('comp1').mesh('mesh1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size2').set('hmin', 6.0E-4);
model.component('comp1').mesh('mesh1').feature('size2').set('hminactive', false);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');

model.component('comp1').mesh('mesh1').create('bl1', 'BndLayer');
model.component('comp1').mesh('mesh1').feature('bl1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('bl1').selection.set([1 2]);
model.component('comp1').mesh('mesh1').feature('bl1').create('blp', 'BndLayerProp');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').selection.set([8 9]);
model.component('comp1').mesh('mesh1').feature('bl1').label('Boundary Layers interface');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').set('inittype', 'blhmin');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').set('blhmin', '0.0002');

model.component('comp1').mesh('mesh1').run;

% study
model.study.create('std2');
model.study('std2').create('stat', 'Stationary');
model.study('std2').feature('stat').set('activate', {'spf2' 'off' 'dl2' 'off' 'tds2' 'on' 'frame:spatial2' 'on' 'frame:material2' 'on'});
model.sol.create('sol4');
model.sol('sol4').study('std2');
model.sol('sol4').attach('std2');
model.sol('sol4').create('st1', 'StudyStep');
model.sol('sol4').create('v1', 'Variables');
model.sol('sol4').create('s1', 'Stationary');
model.sol('sol4').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol4').feature('s1').create('d1', 'Direct');
model.sol('sol4').feature('s1').create('i1', 'Iterative');
model.sol('sol4').feature('s1').create('i2', 'Iterative');
model.sol('sol4').feature('s1').create('i3', 'Iterative');
model.sol('sol4').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol4').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol4').feature('s1').feature('i3').create('mg1', 'Multigrid');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').create('sc1', 'SCGS');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').create('sc1', 'SCGS');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol4').feature('s1').feature.remove('fcDef');

model.study('std2').feature('stat').label('Stationary NS');
model.study('std2').feature('stat').set('usesol', true);
model.study('std2').feature('stat').set('notsolmethod', 'sol');

model.sol('sol4').attach('std2');
model.sol('sol4').feature('st1').label('Compile Equations: Stationary NS');
model.sol('sol4').feature('v1').label('Dependent Variables 1.1');
model.sol('sol4').feature('v1').set('notsolmethod', 'sol');
model.sol('sol4').feature('s1').label('Stationary Solver 1.1');
model.sol('sol4').feature('s1').feature('dDef').label('Direct 2');
model.sol('sol4').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol4').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol4').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol4').feature('s1').feature('fc1').set('linsolver', 'd1');
model.sol('sol4').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol4').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol4').feature('s1').feature('fc1').set('maxiter', 100);
model.sol('sol4').feature('s1').feature('d1').label('Direct, pressure (dl) (merged)');
model.sol('sol4').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol4').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol4').feature('s1').feature('i1').label('AMG, concentrations (tds)');
model.sol('sol4').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol4').feature('s1').feature('i1').set('maxlinit', 1000);
model.sol('sol4').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol4').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol4').feature('s1').feature('i2').label('AMG, pressure (dl)');
model.sol('sol4').feature('s1').feature('i2').set('nlinnormuse', true);
model.sol('sol4').feature('s1').feature('i2').set('maxlinit', 1000);
model.sol('sol4').feature('s1').feature('i2').feature('ilDef').label('Incomplete LU 1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').label('Multigrid 1.1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').set('prefun', 'saamg');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').set('usesmooth', false);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol4').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol4').feature('s1').feature('i3').label('AMG, fluid flow variables (spf)');
model.sol('sol4').feature('s1').feature('i3').set('nlinnormuse', true);
model.sol('sol4').feature('s1').feature('i3').set('maxlinit', 1000);
model.sol('sol4').feature('s1').feature('i3').set('rhob', 20);
model.sol('sol4').feature('s1').feature('i3').feature('ilDef').label('Incomplete LU 1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').label('Multigrid 1.1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').set('prefun', 'saamg');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').set('maxcoarsedof', 80000);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').set('strconn', 0.02);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').set('saamgcompwise', true);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').set('usesmooth', false);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').feature('sc1').label('SCGS 1.1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').feature('sc1').set('approxscgs', true);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('pr').feature('sc1').set('scgsdirectmaxsize', 1000);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').feature('sc1').label('SCGS 1.1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').feature('sc1').set('iter', 1);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').feature('sc1').set('approxscgs', true);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('po').feature('sc1').set('scgsdirectmaxsize', 1000);
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol4').feature('s1').feature('i3').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol4').runAll;

% results
model.result.dataset.create('mesh1', 'Mesh');
model.result.dataset('mesh1').set('mesh', 'mesh1');

model.result.dataset.create('dset2', 'Solution');
model.result.dataset('dset2').label('Solution MCTC 2D');

model.result.dataset.create('edg2', 'Edge2D');
model.result.dataset('edg2').set('data', 'dset2');
model.result.dataset('edg2').selection.set([8 9]);
model.result.dataset('edg2').label('Interface 2D');

model.result.dataset.create('edg3', 'Edge2D');
model.result.dataset('edg3').set('data', 'dset2');
model.result.dataset('edg3').selection.set([4 5]);
model.result.dataset('edg3').label('MCTC wall');

model.result.dataset.create('edg4', 'Edge2D');
model.result.dataset('edg4').set('data', 'dset2');
model.result.dataset('edg4').selection.set([3]);
model.result.dataset('edg4').label('Top Wall');

% model.result.export.create('data1', 'Data');
% model.result.export('data1').label('Interface');
% model.result.export('data1').set('data', 'edg2');
% model.result.export('data1').set('expr', {'c2'});
% model.result.export('data1').set('unit', {'mol/m^3'});
% model.result.export('data1').set('descr', {'Concentration'});
% model.result.export('data1').set('filename', 'C:\Users\Willy\Work\Comsol_simulations\MCTC_reservoir_rigid\2D\50pc_sphere_unconfined\output_interface\out.txt');
% model.result.export('data1').run;

model.result.export.create('data2', 'Data');
model.result.export('data2').label('All');
model.result.export('data2').set('expr', {'c2'});
model.result.export('data2').set('unit', {'mol/m^3'});
model.result.export('data2').set('descr', {'Concentration'});
model.result.export('data2').set('filename', 'C:\Users\Willy\Work\Comsol_simulations\MCTC_reservoir_rigid\2D\50pc_sphere_unconfined\output_isolated\out.txt');
model.result.export('data2').run;

data = load('C:\Users\Willy\Work\Comsol_simulations\MCTC_reservoir_rigid\2D\50pc_sphere_unconfined\output_isolated\out.txt');
x = data(:,1);
y = data(:,2);
c = data(:,3);
save(['C:\Users\Willy\Work\Comsol_simulations\MCTC_reservoir_rigid\2D\50pc_sphere_unconfined\output_isolated\' ...
    'Da_' num2str(argin(3)) '.mat'],'x','y','c');

out = model;

