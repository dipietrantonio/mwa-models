#!/usr/bin/env python3
from PIL import Image
import sys


def draw_rect(pixels, at, width, height):
    selected_coords = []
    start_row, start_col = at
    line_width = 2
    for row in range(start_row, start_row + height + line_width):
        for w in range(line_width):
            selected_coords.append((row, start_col + w))
            selected_coords.append((row, start_col + width + w))
    
    for col in range(start_col, start_col + width):
        for w in range(line_width):
            selected_coords.append((start_row + w, col))
            selected_coords.append((start_row + w + height, col))
    
    for i, j in selected_coords:
        pixels[i,j] = (255, 255, 255, 255)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    im = Image.open(input_file)
    width, height = im.size
    pixels = im.load() # create the pixel map
    side = 20
    draw_rect(pixels, (85, 115), side, side)
    im.save(output_file, "png", quality=100)
    # im_cropped = im.crop(margins)
    # im_cropped = im_cropped.resize((1024, 1024))
    # im_cropped.show()
    # im_cropped.save(output_file, "png", quality=95)


