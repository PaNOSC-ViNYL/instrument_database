project('mcpl',
        'c',
        version: '1.6.2',
        #        #licence: 'GPL v2',
        default_options: ['default_library=static']
       )

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
zlib_dep = dependency('zlib')
public_headers_inc = include_directories('src/mcpl')
mcpl_lib  = library('mcpl',
                    version: meson.project_version(),
                    soversion: '0',
                    #cpp_args: ['-Wno-comment'],
                    sources: [ 'src/mcpl/mcpl.c' ],
                    include_directories: public_headers_inc,
                    dependencies: [m_dep,zlib_dep],
                    #install: true,
                    #install_tag: 'client',
                   )

############################### Dependencies
mcpl_dep = declare_dependency(
  #sources: cameo_main_header_file,
  link_with: mcpl_lib,
  include_directories: [public_headers_inc, meson.current_build_dir()],
  #dependencies: m_dep,
)

if not meson.is_subproject()
############################### pkgconfig
pkg_mod = import('pkgconfig')
pkg_mod.generate(mcpl_lib,
                 version: meson.project_version(),
                 #name: 'libcameo-api-cpp',
                 #filebase: 'cameo-api-cpp',
                 description: 'MCPL library',
#                 subdirs: include_path_subdir,
                )

install_headers('src/mcpl/mcpl.h')
endif
