project('dilutedlat', ['cpp', 'c'],
  version: '1.0.0',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=debug', 'cpp_std=c++2a' ]
)

if get_option('buildtype') == 'release'
  add_project_arguments('-Ofast',  language : 'cpp')
endif

latlib_dep = dependency('liblatindex')

diluter_bin = executable('dmnd_dilute', 
  files('dmnd_dilute.cpp'),
  dependencies: latlib_dep,
  include_directories: 'include'
  )

#diluter_nd_bin = executable('dmnd_dilute_nodelete',
#  files('dmnd_dilute_nodelete.cpp'),
#  dependencies: latlib_dep,
#  include_directories: 'include'
#  )
