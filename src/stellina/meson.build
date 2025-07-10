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

# Use the same dependencies as the main Siril build
stellina_deps = [
  cc.find_library('m', required: false),
  dependency('json-glib-1.0', required: true),
  dependency('glib-2.0', required: true),
  dependency('gtk+-3.0', required: true),  # Add GTK dependency
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

# Export the dependency for use in main siril executable
stellina_dep = declare_dependency(
  link_with: libstellina,
  include_directories: include_directories('.'),
  dependencies: stellina_deps
)