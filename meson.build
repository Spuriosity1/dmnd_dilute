project('dilutedlat', ['cpp', 'c'],
  version: '1.0.0',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=release', 'cpp_std=c++2a' ]
)

if get_option('buildtype') == 'release'
  add_project_arguments('-Ofast',  language : 'cpp')
endif

# Prefers a system install (via pkg-config); otherwise builds latticelab as a
# meson subproject from subprojects/latticelab.wrap
latlib_dep = dependency('liblatindex', version : '>=1.1',
  required: true,
  fallback : ['latticelab', 'latlib_dep'],
  default_options : ['enable-tests=false', 'enable-benchmarks=false'],
)
json_dep = dependency('nlohmann_json', required: true)

diluter_bin = executable('dmnd_dilute', 
  files('src/dmnd_dilute.cpp'),
  dependencies: [latlib_dep,
      json_dep
    ],
  include_directories: 'include'
  )

count_neighbours_bin = executable('count_neighbours',
  files('src/count_neighbours.cpp'),
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
