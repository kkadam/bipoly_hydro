  100 format('********************************',/,                     &
             'E  model_type = ',i2,'       E',/,                       &
             'E  not supported in this       E',/,                     &
             'E  version                     E',/,                     &
             'E         STOPPING             E',/,                     &
             '********************************',///) 
 
  110 format('>>>> This simulation starts from Initial  <<<<',/,       &
             '>>>> Conditions from an SCF initial model <<<<',///)

  120 format('>>>> This simulation starts from Initial  <<<<',/,       &
             '>>>> Conditions from a continuation file  <<<<',///)

  130 format('************************************************',/,     &
             '*     This simulation assumes No Symmetries    *',/,     &
             '*                                              *',/,     &
             '* Evolution will proceed from timestep: ',i7,'*',/,      &
             '*                                              *',/,     &
             '*                          to timestep: ',i7,'*',/,      &
             '************************************************',///)

  140 format('************************************************',/,     &
             '* This simulation assumes Equatorial Symmetrty *',/,     &
             '*                                              *',/,     &
             '* Evolution will proceed from timestep: ',i6,'*',/,      &
             '*                                              *',/,     &
             '*                          to timestep: ',i6,'*',/,      &
             '************************************************',///)

  150 format('************************************************',/,     &
             '*     This simulation assumes Pi Symmetry      *',/,     &
             '*                                              *',/,     &
             '* Evolution will proceed from timestep: ',i6,'*',/,      &
             '*                                              *',/,     &
             '*                          to timestep: ',i6,'*',/,      &
             '************************************************',///)

  160 format('Radial resolution: ',i4,' zones',/,                      &
             'Vertical resolution: ',i4,' zones',/,                    &
             'Azimuthal resolution: ',i4,' zones',//)        

  170 format('Polytropic Index',T20,'Polytropic Exponent',T40,         &
             'Kappa *1',T60,'Kappa *2',/,                              &
             f12.5,T20,f12.5,T40,f12.5,T60,f12.5,//)

  180 format('Vacuum Density',T20,'Vacuum Tau',T40,'Background Pressure',/,    &
             es12.5,T20,es12.5,T40,f12.5,//)

  190 format('Timstep will be set to ',f6.3,' times the smallest',/,           &
             'signal crossing time',//)

  200 format('Velocities of material less dense than ',es12.5,/,               &
             'will be limited to ',f6.3,' times the maximum sound speed',/,    &
             'Velocities of matrerial more dense than the cutoff will',/,      &
             'be limited to ',es12.5,//)

  210 format(f12.5,' frames will be written per orbit',//,             &
             'the orbital period is ',f12.5,//,                        &
             'Diagnostics written every ',i4,' timestep',//)

  220 format('boundary density for definition of a star is ',f12.5,//,         &
             'boundary point lies at ',f12.5,' from system center of mass to density max of primary',//)

  225 format('coefficient of artificial viscosity is ',f12.5,//)

  230 format('>>>> Simulation Includes Forces from Self-Gravity of Fluid <<<<',//)


