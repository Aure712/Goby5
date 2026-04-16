from PIL import Image
Image.MAX_IMAGE_PIXELS = None
# === 配置区域 ============================================
# 输入图片路径（改成你的图片）
input_image_path = "Tra-1b_RT.png"

# 输出文件名（可自行修改）
output_crop1_path = "HE_region_1.png"
output_crop2_path = "HE_region_2.png"

# 裁剪区域 1（左上角 6762, 15111，宽 4830，高 2898）
x1, y1 = 6762, 15111
w1, h1 = 4830, 2898
crop_box_1 = (x1, y1, x1 + w1, y1 + h1)

# 裁剪区域 2（左上角 9062, 5451，宽 4830，高 2898）
x2, y2 = 9062, 5451
w2, h2 = 4830, 3611
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
