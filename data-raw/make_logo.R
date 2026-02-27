# ============================================================================
# make_logo.R --- Generate hurdlebb hex sticker logo (PNG from SVG)
#
# The logo is a hand-crafted SVG at man/figures/logo.svg.
# This script converts the SVG to a 500px-wide PNG for use in README
# and pkgdown. Edit the SVG directly for design changes.
#
# Requirements: magick (with rsvg support)
# ============================================================================

library(magick)

svg_path <- file.path("man", "figures", "logo.svg")
png_path <- file.path("man", "figures", "logo.png")

svg <- image_read_svg(svg_path, width = 500)
image_write(svg, png_path, format = "png")

cat("Logo PNG generated:", png_path, "\n")
cat("Size:", file.info(png_path)$size, "bytes\n")
