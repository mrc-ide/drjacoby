# Create package hex image.
path_image <- "R_ignore/images/Sunglasses.png"
s <- hexSticker::sticker(subplot = path_image,
                         package = "Dr Jacoby", p_size = 80, p_color = "darkred",
                         s_x = 1, s_y = 0.8, s_width = 0.8, s_height = 0.8,
                         h_fill = "white",
                         h_color = "black",
                         dpi = 1000,
                         filename = "R_ignore/images/logo3.png",
                         white_around_sticker = F)
