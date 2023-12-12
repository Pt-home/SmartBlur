#include <emscripten.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

typedef unsigned char byte;
typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef float float32;

void applyBoxBlur(byte* image, int width, int height, int radius) {
    int size = radius * 2 + 1; // Размер окна размытия

    byte* tempData = (byte*)malloc(width * height * 4);
    memcpy(tempData, image, width * height * 4);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int sumR = 0, sumG = 0, sumB = 0;

            // Проходимся по окну размытия
            for (int i = -radius; i <= radius; i++) {
                for (int j = -radius; j <= radius; j++) {
                    int neighborX = x + j;
                    int neighborY = y + i;

                    if (neighborX >= 0 && neighborX < width && neighborY >= 0 && neighborY < height) {
                        int neighborIndex = (neighborY * width + neighborX) * 4;
                        sumR += tempData[neighborIndex];
                        sumG += tempData[neighborIndex + 1];
                        sumB += tempData[neighborIndex + 2];
                    }
                }
            }

            // Вычисляем средние значения цветов
            int avgR = sumR / (size * size);
            int avgG = sumG / (size * size);
            int avgB = sumB / (size * size);

            // Обновляем цвет пикселя
            int index = (y * width + x) * 4;
            image[index] = avgR;
            image[index + 1] = avgG;
            image[index + 2] = avgB;
        }
    }

    free(tempData);
}

// Helper function to apply Gaussian blur to the cloned image
void applyGaussianBlur(byte* image, int width, int height, int radius) {
    // Apply Gaussian blur operation to the image using your preferred method
    // You can use algorithms like Gaussian blur or separable Gaussian blur
    // Here's a simplified example using a separable Gaussian blur for demonstration purposes:

    byte* temp = (byte*)malloc(height * width * 4);
    memcpy(temp, image, height * width * 4);

    float* kernel = (float*)malloc((2 * radius + 1) * sizeof(float));

    // Generate Gaussian kernel
    float sigma = radius / 3.0f;
    float twoSigmaSquare = 2.0f * sigma * sigma;
    float sqrtTwoPiSigmaSquare = sqrtf(2.0f * M_PI * sigma * sigma);
    float sum = 0.0f;

    for (int i = -radius; i <= radius; i++) {
        float exponent = -((float)(i * i)) / twoSigmaSquare;
        kernel[i + radius] = expf(exponent) / sqrtTwoPiSigmaSquare;
        sum += kernel[i + radius];
    }

    // Normalize the kernel
    for (int i = 0; i < 2 * radius + 1; i++) {
        kernel[i] /= sum;
    }

    // Apply horizontal blur
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int pixelIndex = (y * width + x) * 4;

            float sumR = 0.0f, sumG = 0.0f, sumB = 0.0f;

            for (int i = -radius; i <= radius; i++) {
                int newX = x + i;

                if (newX >= 0 && newX < width) {
                    int newIndex = (y * width + newX) * 4;
                    float weight = kernel[i + radius];
                    sumR += weight * temp[newIndex];
                    sumG += weight * temp[newIndex + 1];
                    sumB += weight * temp[newIndex + 2];
                }
            }

            image[pixelIndex] = (byte)sumR;
            image[pixelIndex + 1] = (byte)sumG;
            image[pixelIndex + 2] = (byte)sumB;
        }
    }

    // Apply vertical blur
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int pixelIndex = (y * width + x) * 4;

            float sumR = 0.0f, sumG = 0.0f, sumB = 0.0f;

            for (int i = -radius; i <= radius; i++) {
                int newY = y + i;

                if (newY >= 0 && newY < height) {
                    int newIndex = (newY * width + x) * 4;
                    float weight = kernel[i + radius];
                    sumR += weight * image[newIndex];
                    sumG += weight * image[newIndex + 1];
                    sumB += weight * image[newIndex + 2];
                }
            }

            temp[pixelIndex] = (byte)sumR;
            temp[pixelIndex + 1] = (byte)sumG;
            temp[pixelIndex + 2] = (byte)sumB;
        }
    }

    // Copy the blurred image back to the original image buffer
    memcpy(image, temp, height * width * 4);

    // Clean up temporary buffers
    free(temp);
    free(kernel);
}

// Helper function to apply surface blur to the image
void applySurfaceBlur(byte* image, int width, int height, int radius) {
    // Apply surface blur operation to the image using your preferred method
    // You can use algorithms like the mean shift or anisotropic diffusion
    // Here's a simplified example using a mean shift algorithm for demonstration purposes:

    byte* temp = (byte*)malloc(height * width * 4);
    memcpy(temp, image, height * width * 4);

    // Define the kernel size based on the radius
    int kernelSize = radius * 2 + 1;

    // Apply mean shift filter to each pixel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int pixelIndex = (y * width + x) * 4;

            int sumR = 0, sumG = 0, sumB = 0, count = 0;

            // Compute the mean shift for the pixel
            for (int i = -radius; i <= radius; i++) {
                for (int j = -radius; j <= radius; j++) {
                    int newY = y + i;
                    int newX = x + j;

                    if (newY >= 0 && newY < height && newX >= 0 && newX < width) {
                        int newIndex = (newY * width + newX) * 4;
                        sumR += temp[newIndex];
                        sumG += temp[newIndex + 1];
                        sumB += temp[newIndex + 2];
                        count++;
                    }
                }
            }

            // Compute the average color
            byte avgR = (byte)(sumR / count);
            byte avgG = (byte)(sumG / count);
            byte avgB = (byte)(sumB / count);

            // Update the pixel in the image buffer
            image[pixelIndex] = avgR;
            image[pixelIndex + 1] = avgG;
            image[pixelIndex + 2] = avgB;
        }
    }

    // Clean up temporary buffer
    free(temp);
}

void findEdges(unsigned char* imageData, unsigned char* edgeMask, int width, int height, float threshold) {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = (y * width + x) * 4;

            unsigned char r = imageData[index];
            unsigned char g = imageData[index + 1];
            unsigned char b = imageData[index + 2];

            // Вычисляем сумму разниц цветов текущего пикселя с его соседями
            float colorDiffSum = 0.0;
            int neighborCount = 0;

            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    int neighborX = x + i;
                    int neighborY = y + j;

                    if (neighborX >= 0 && neighborX < width && neighborY >= 0 && neighborY < height &&
                        !(i == 0 && j == 0)) {
                        int neighborIndex = (neighborY * width + neighborX) * 4;
                        unsigned char neighborR = imageData[neighborIndex];
                        unsigned char neighborG = imageData[neighborIndex + 1];
                        unsigned char neighborB = imageData[neighborIndex + 2];

                        float colorDiff = sqrtf(powf(r - neighborR, 2) + powf(g - neighborG, 2) + powf(b - neighborB, 2));
                        colorDiffSum += colorDiff;
                        neighborCount++;
                    }
                }
            }

            // Вычисляем среднее значение разницы цветов
            float avgColorDiff = colorDiffSum / neighborCount;

            // Пороговое значение для определения краев
            if (avgColorDiff > threshold) {
                edgeMask[index] = 255; // Устанавливаем белый цвет (край)
                edgeMask[index + 1] = 255;
                edgeMask[index + 2] = 255;
            } else {
                edgeMask[index] = 0; // Устанавливаем черный цвет (не край)
                edgeMask[index + 1] = 0;
                edgeMask[index + 2] = 0;
            }
        }
    }
}

float calculateDistanceToEdge(unsigned char* maskData, int width, int height, int x, int y, int size) {
    float minDistance = width * height;
    int center = size / 2;

    for (int i = -center; i <= center; i++) {
        for (int j = -center; j <= center; j++) {
            int neighborX = x + i;
            int neighborY = y + j;

            if (neighborX >= 0 && neighborX < width && neighborY >= 0 && neighborY < height) {
                int neighborIndex = (neighborY * width + neighborX) * 4;
                if (maskData[neighborIndex] == 255) {
                    float distance = sqrtf((i * i) + (j * j));
                    if (distance < minDistance) {
                        minDistance = distance;
                    }
                }
            }
        }
    }

    return (float)minDistance;
}

void interpolateColors(unsigned char* imageData, unsigned char* originalData, unsigned char* maskData, int width, int height, int size, const char* mode) {
    // Вычисляем максимальное возможное расстояние в окрестности
    float maxDistance = sqrtf(2 * size * size);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = (y * width + x) * 4;

            unsigned char originalR = originalData[index];
            unsigned char originalG = originalData[index + 1];
            unsigned char originalB = originalData[index + 2];
            unsigned char originalA = originalData[index + 3];

            unsigned char blurredR = imageData[index];
            unsigned char blurredG = imageData[index + 1];
            unsigned char blurredB = imageData[index + 2];
            unsigned char blurredA = imageData[index + 3];

            // Проверяем, является ли пиксель краем
            if (maskData[index] == 255) {
				if (strcmp(mode, "normal") == 0) {
					// Пиксель на краю - оставляем его без изменений
					imageData[index]     = originalR;
					imageData[index + 1] = originalG;
					imageData[index + 2] = originalB;
					imageData[index + 3] = originalA;
				} else {
					// Пиксель на краю - overlay
					imageData[index]     = 255;
					imageData[index + 1] = 255;
					imageData[index + 2] = 255;
					imageData[index + 3] = 255;
				}
                
            } else if (size > 0) {
                // Вычисляем расстояние до края (относительное расстояние)
                float distance = calculateDistanceToEdge(maskData, width, height, x, y, size);
				if (distance < width*height-1) {
					// Вычисляем прозрачность размытому
					float ratio = distance / maxDistance;
					unsigned char distancedA = roundf(ratio * 255);
					
					//Нормируем
					float r1_norm = blurredR / 255.0f; 
					float g1_norm = blurredG / 255.0f;
					float b1_norm = blurredB / 255.0f;
					float a1_norm = distancedA / 255.0f;

					float r2_norm = originalR / 255.0f;
					float g2_norm = originalG / 255.0f;
					float b2_norm = originalB / 255.0f;
					float a2_norm = originalA / 255.0f;

					//Итоговый цвет (в нормализованном формате):
					float r_result_norm = a1_norm * r1_norm + (1 - a1_norm) * a2_norm * r2_norm;
					float g_result_norm = a1_norm * g1_norm + (1 - a1_norm) * a2_norm * g2_norm;
					float b_result_norm = a1_norm * b1_norm + (1 - a1_norm) * a2_norm * b2_norm;
					float a_result_norm = a1_norm + (1 - a1_norm) * a2_norm;

					// Обновляем цвет и прозрачность пикселя
					imageData[index]     = roundf(r_result_norm * 255);
					imageData[index + 1] = roundf(g_result_norm * 255);
					imageData[index + 2] = roundf(b_result_norm * 255);
					imageData[index + 3] = roundf(a_result_norm * 255);
				}
            }
        }
    }
}



// Главная функция применения Smart Blur с различными методами
void applySmartBlurWithType(unsigned char* imageData, int width, int height, float radius, float threshold, int size, const char* type, const char* mode) {
    unsigned char* originalData = (unsigned char*)malloc(width * height * 4);
    memcpy(originalData, imageData, width * height * 4);

    // Применяем размытие (в зависимости от типа)
    if (strcmp(type, "gaussian") == 0) {
        applyGaussianBlur(imageData, width, height, radius);
    } else if (strcmp(type, "surface") == 0) {
        applySurfaceBlur(imageData, width, height, radius);
    } else if (strcmp(type, "box") == 0) {
        applyBoxBlur(imageData, width, height, radius);
	} 
	
    // Создаем маску краев
    unsigned char* maskData = (unsigned char*)malloc(width * height * 4);
    findEdges(originalData, maskData, width, height, threshold);

    // Интерполируем цвета на основе маски краев
    interpolateColors(imageData, originalData, maskData, width, height, size, mode);

    free(originalData);
    free(maskData);
}

// Главная функция применения Smart Blur с различными методами
void applySmartBlurEdge(unsigned char* imageData, int width, int height, float threshold) {
    
    // Создаем маску краев
    unsigned char* maskData = (unsigned char*)malloc(width * height * 4);
    findEdges(imageData, maskData, width, height, threshold);
	memcpy(imageData, maskData, width * height * 4);
    free(maskData);
	
}


EMSCRIPTEN_KEEPALIVE
byte* wasmAlloc(uint32 width, uint32 height) {
  return malloc(width * height * 4);
}

EMSCRIPTEN_KEEPALIVE
void wasmFree(byte* p) {
  free(p);
}

EMSCRIPTEN_KEEPALIVE
void wasmProcess(byte* imageData, int width, int height, float radius, float threshold, const char* type, int size, const char* mode) {
	
  if (strcmp(mode, "edge") == 0) {        
		applySmartBlurEdge(imageData, width, height, threshold);
    } else {
        applySmartBlurWithType(imageData, width, height, radius, threshold, size, type, mode);
    }
  
}

