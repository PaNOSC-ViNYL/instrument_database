project('instrument_repo', 'c',
        license: 'EUPL-1.2',
       )

py_to_instr = files('mcstas/scripts/py_to_instr.py')
mcstas_custom_components_dir = 'mcstas/components/'
project_source_root = meson.project_source_root()

# configurations
instruments_bindir = 'mcstas_instruments'

instrfiledir = 'instruments'
instruments_dict = { 'D11' : 'simple', 'Panther': 'quick', 'ThALES': 'quick'}
sample_list = ['None', 'H2O', 'D2O', 'vanadium', 'sqw','qSq']
sample_holder_list = ['None', 'QUARTZ']
sample_environment_list = ['None']


python = find_program('python3')
mcstas = find_program('mcstas')
nvc = find_program('nvc', required: false)
if(nvc.found())
  mcstas_gpu_flags = [ '-fast', '-acc=gpu', '-gpu=managed', '-DOPENACC','-Wno-unused', '--diag_suppress=incompatible_param','--diag_suppress=used_before_set']
  openacc_link_flags = ['-lm', '-acc=gpu', '-L/usr/lib/x86_64-linux-gnu']
  gpu_suffix='-gpu'
endif

cc = meson.get_compiler('c')
message(meson.get_compiler('c').get_id())
m_dep = dependency('', required : false)
cpu_suffix=''
mcstas_c_flags = ['-O3', '-Wno-all', '-Wformat-truncation=0', '-Wformat-overflow=0', '-Wno-format-extra-args', '-Wno-format','-Wno-unknown-pragmas']
m_dep = cc.find_library('m', required : false)

mpi_dep = dependency('mpi',language:'c')
if(mpi_dep.found())
  #mcstas_c_flags +=['-DUSE_MPI']
endif

#mcpl_config = find_program('mcpl-config', version: '1.6.2', required: false)
#if(not mcpl_config.found() or true)
mcpl_dep = dependency('mcpl', version: '1.6.2', fallback: ['mcpl'], required: true, static:true)
#endif

DISABLE_GZIP_MCPL = '-DNOGZIP'

sources_list = []
                                      
foreach instrname, flavour : instruments_dict
  foreach sample : sample_list
    foreach sample_holder : sample_holder_list
      foreach sample_environment : sample_environment_list
        
        c = run_command([python, py_to_instr, '--dry_run', 'ILL', instrname, 'HEAD', 'mcstas', flavour, sample, sample_environment, sample_holder], check: true)
        # get the name of the executables for the different stages as a list of strings
        ss = c.stdout().strip().split('\n')[-1].strip('[]').strip().replace('\'','').replace(' ','').strip().split(',')
        #message(ss)
        foreach  s : ss
          if not sources_list.contains(s) # list of already defined sources, to avoid duplicates
            sources_list += s
            outname =s.strip(' \'').strip('\'')
            output = outname.split('/')[-1]
            
            message('Creating: '+output)
            instr = custom_target(output,
                                  build_by_default: true,
                                  command: [python, py_to_instr, '--outdir', '@OUTDIR@',
                                            'ILL', instrname, 'HEAD', 'mcstas', flavour, sample, sample_environment, sample_holder],
                                  install: false,
                                  output: output
                                 )
            output=output.replace('.instr','.c')
            message('Creating: '+output)
            mytarget = custom_target(output,
                                     build_by_default: true,
                                     output: output,
                                     input: instr,
                                     command: [mcstas, '--verbose', '-t', '-I', project_source_root/mcstas_custom_components_dir, '-o', '@OUTPUT@', '@INPUT@'],
                                     install: false
                                    )
            
            exe = executable(output.replace('.c',cpu_suffix+'.out'), mytarget,
                             dependencies : [m_dep, mcpl_dep,mpi_dep],
                             c_args: [ mcstas_c_flags, DISABLE_GZIP_MCPL],
                             #link_args: openacc_link_flags,
                             install_dir: get_option('bindir')+'/'+instruments_bindir,
                             install: true,
                             install_tag: 'instruments',
                             #install_rpath: '$ORIGIN/../../lib' #,/usr/lib/x86_64-linux-gnu'
                            )
            
            if(nvc.found())
              exegpu = custom_target(output.replace('.c',gpu_suffix+'.out'),
                                     input: mytarget,
                                     output: output.replace('.c',gpu_suffix+'.out'),
                                     command: [ 'nvc', mcstas_gpu_flags, openacc_link_flags, '-o', '@OUTPUT@', '@INPUT@', '-L' + mcpl_lib_dir, '-lmcpl', '-Wl,-rpath,'+mcpl_lib_dir],
                                     depends : mytarget,
                                     install: true,
                                     install_tag: 'instruments',
                                     install_dir: get_option('bindir')+'/'+instruments_bindir,
                                     #install_rpath: ',/usr/lib/x86_64-linux-gnu',
                                    )
              
              xmlfile = configure_file(configuration: {'exe': output.replace('.c',gpu_suffix),
                                                       'instruments_bindir': instruments_bindir,
                                                       'bindir': get_option('prefix')+'/'+get_option('bindir'),
                                                       'exesuffix': '',
                                                      },
                                       input: files('xml_file.in'),
                                       output: output.replace('.c', gpu_suffix+'.xml'),
                                       install: false
                                      )
              xmlfiles+=xmlfile
            endif
            
          endif
          
        endforeach
      endforeach
    endforeach
  endforeach
endforeach

