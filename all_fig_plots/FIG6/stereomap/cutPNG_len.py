from PIL import Image
Image.MAX_IMAGE_PIXELS = None
# === 配置区域 ============================================
# 输入图片路径（改成你的图片）
input_image_path = "FUN_017191_-3844_18544_0_11800.png"

# 输出文件名（可自行修改）
output_crop1_path = "FUN_017191_count.png"
output_crop2_path = "FUN_017191_bar.png"

# 裁剪区域 1（左上角 6875, 5126，宽 3146，高 1452）
x1, y1 = 15488, 7744
w1, h1 = 440, 825
crop_box_1 = (x1, y1, x1 + w1, y1 + h1)

# 裁剪区域 2（左上角 8404, 682，宽 3146，高 2068）
x2, y2 = 3364, 8244
w2, h2 = 766, 312
crop_box_2 = (x2, y2, x2 + w2, y2 + h2)
# ==========================================================


def crop_image(img_path, box, save_path):
    img = Image.open(img_path)
    cropped = img.crop(box)
    cropped.save(save_path)
    print(f"已输出：{save_path}")


if __name__ == "__main__":
    # 裁剪区域 1
    crop_image(input_image_path, crop_box_1, output_crop1_path)
    
    # 裁剪区域 2
    crop_image(input_image_path, crop_box_2, output_crop2_path)

    print("全部裁剪完成！")
