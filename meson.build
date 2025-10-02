project('dilutedlat', ['cpp', 'c'],
  version: '1.0.0',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=release', 'cpp_std=c++2a' ]
)

if get_option('buildtype') == 'release'
  add_project_arguments('-Ofast',  language : 'cpp')
endif

latlib_dep = dependency('liblatindex', version : '>=1.1', required: true)
json_dep = dependency('nlohmann_json', required: true)

#latlib_proj = subproject('liblatindex')
#latlib_dep = latlib_proj.get_variable('latlib_dep')

diluter_bin = executable('dmnd_dilute', 
  files('dmnd_dilute.cpp'),
  dependencies: [latlib_dep,
      json_dep
    ],
  include_directories: 'include'
  )

#diluter_nd_bin = executable('dmnd_dilute_nodelete',
#  files('dmnd_dilute_nodelete.cpp'),
#  dependencies: latlib_dep,
#  include_directories: 'include'
#  )
