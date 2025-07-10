# src/stellina/meson.build
# Stellina processing extension for Siril

stellina_sources = [
  'stellina_processor.cpp',
  'stellina_commands.cpp',
  'stellina_coordinates.cpp',
  'stellina_quality.cpp'
]

# Direct libnova detection - we know it's installed at /opt/homebrew
libnova_inc = include_directories('/opt/homebrew/include')
libnova_lib = cc.find_library('nova', dirs: ['/opt/homebrew/lib'], required: true)

libnova_dep = declare_dependency(
  include_directories: libnova_inc,
  dependencies: libnova_lib
)

stellina_deps = [
  cc.find_library('m', required: false),
  dependency('json-glib-1.0', required: true),
  libnova_dep,
  cfitsio_dep,
  gsl_dep
]

# Create stellina static library
libstellina = static_library(
  'stellina',
  stellina_sources,
  dependencies: stellina_deps,
  include_directories: [siril_inc_dir, include_directories('.')],
  install: false
)

# Add to main Siril dependencies
stellina_dep = declare_dependency(
  link_with: libstellina,
  include_directories: include_directories('.'),
  dependencies: stellina_deps
)