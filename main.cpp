#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <fstream>

// if you use linux, windows 10, windows 11 enable use_ansi this flag, in other case (windows 8, 7... set to 0)
#define use_ansi 1

typedef struct Size { int width, height; } Size;
typedef struct vec3 { double x, y, z; } vec3;
typedef enum ColorBit { Black, DarkBlue, DarkGreen, DarkCyan, DarkRed, DarkMagenta, DarkYellow, Gray, DarkGray, Blue, Green, Cyan, Red, Magenta, Yellow, White } ColorBit;

inline double fast_atan2(double y, double x) {
	if (x == 0.0) {
		if (y > 0.0) { return 1.57079632679; }
		else if (y < 0.0) { return -1.57079632679; }
		return 0.0;
	}
	double angle, xx, div = y / x;
	if (div < -1.0 || div > 1.0) {
		div = 1.0 / div;
		xx = div * div;
		angle = ((0.0776509570923569 * xx + -0.287434475393028) * xx + 0.995181681698) * div;
		angle = div > 0.0 ? (1.57079632679 - angle) : (-1.57079632679 - angle);
		if (x < 0.0) { angle += (y >= 0.0 ? 3.141592653589793 : -3.141592653589793); }
	}
	else {
		xx = div * div;
		angle = ((0.0776509570923569 * xx + -0.287434475393028) * xx + 0.995181681698) * div;
		if (x < 0.0) { angle += (y >= 0.0 ? 3.141592653589793 : -3.141592653589793); }
	}
	return angle;
}

inline double fast_cos(double x) {
	x *= 0.159154943092;
	x -= 0.25 + floor(x + 0.25);
	x *= 16.0 * (fabs(x) - 0.5);
	return 0.224 * x * (fabs(x) - 1.0) + x;
}

inline double fast_sin(double x) {
	x = (x - 1.57079632679) * 0.159154943092;
	x -= 0.25 + floor(x + 0.25);
	x *= 16.0 * (fabs(x) - 0.5);
	return 0.224 * x * (fabs(x) - 1.0) + x;
}

inline double fast_sqrt(double x) {
	double x2 = x * 0.5, y = x;
	long long i = *(long long*)&y;
	i = 0x5fe6ec85e7de30da - (i >> 1);
	y = *(double*)&i;
	return x * y * (1.5 - (x2 * y * y));
}

inline double fast_pow(double a, double b) {
	if (b > 0.0) {
		if (b == 2.0) return a * a;
		if (b == 0.5) return fast_sqrt(a);
		if (b == 1.5) return fast_sqrt(a) * a;
		return pow(a, b);
	}
	if (b == -1.0) return 1.0 / a;
	if (b == -0.5) return 1.0 / fast_sqrt(a);
	if (b == -1.5) return 1.0 / (fast_sqrt(a) * a);
	return 1.0 / pow(a, -b);
}

inline double fast_log(double a) {
	union { double d; long long i; } u;
	u.d = a;
	return (u.i - 4606921278410026770LL) * 1.539095918623324e-16;
}

inline int clamp(int x, int minimum, int maximum) {
	return x <= minimum ? minimum : x >= maximum ? maximum : x;
}

inline int bounce(int x, int min, int max, int divide_max_value) {
	if (x >= min && x <= max) { return x; }
	int offset = (x / divide_max_value) % max;
	return (x / max % 2) ? max - offset : offset;
}

#define min(a, b) a < b ? a : b
#define atan2 fast_atan2
#define cos fast_cos
#define sin fast_sin

#ifdef _WIN32
#include <Windows.h>
#include <conio.h>

HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
DWORD written = 0;

inline void gotoxy(int x, int y) {
	COORD coord = { (short)x, (short)y };
	SetConsoleCursorPosition(h, coord);
}

inline void clear_screen() {
	system("cls");
}

inline Size get_size_screen() {
	CONSOLE_SCREEN_BUFFER_INFO screen;
	GetConsoleScreenBufferInfo(h, &screen);
	return Size { (screen.srWindow.Right - screen.srWindow.Left + 1), (screen.srWindow.Bottom - screen.srWindow.Top + 1) };
}

inline void color(ColorBit font, ColorBit bg) {
	SetConsoleTextAttribute(h, ((bg % 16) * 16) + (font % 16));
}

inline int __kbhit() {
	return _kbhit();
}

inline int getch_toupper() {
	int button = toupper(_getch());
	if (button == 0 || button == 224) {
		unsigned char new_button = _getch();
		if (new_button == 72) { button = 'W'; }
		if (new_button == 80) { button = 'S'; }
		if (new_button == 75) { button = 'A'; }
		if (new_button == 77) { button = 'D'; }
	}
	return button;
}

inline void print_console(const char* str, int length) {
	WriteConsoleA(h, str, length, &written, 0);
}
#endif // Windows

#ifdef __linux__
#include <stdio.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <map>

inline void gotoxy(int x, int y) {
	std::cout << "\033[" << (y + 1) << "d\033[" << (x + 1) << "G" << std::flush;
}

inline void clear_screen() {
	std::cout << "\033c" << std::flush;
}

inline Size get_size_screen() {
	struct winsize screen;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &screen);
	return Size { screen.ws_col, screen.ws_row };
}

inline void color(ColorBit font, ColorBit bg) {
	static const int ansi_to_windows_color[16] = { 0, 4, 2, 6, 1, 5, 3, 7, 8, 12, 10, 14, 9, 13, 11, 15 };
	std::cout << "\033[38;5;" << ansi_to_windows_color[(int)font] << "m\033[48;5;" << ansi_to_windows_color[(int)bg] << "m";
}

int __getch() {
	struct termios old_terminal, new_terminal;
	tcgetattr(STDIN_FILENO, &old_terminal);
	new_terminal = old_terminal;
	new_terminal.c_lflag &= ~(ICANON | ECHO);
	tcsetattr(STDIN_FILENO, TCSANOW, &new_terminal);
	int ch = getchar();
	tcsetattr(STDIN_FILENO, TCSANOW, &old_terminal);
	return ch;
}

int __kbhit() {
	static const int STDIN = 0;
	static int initialized = 0;
	if (!initialized) {
		termios term;
		tcgetattr(STDIN, &term);
		term.c_lflag &= ~ICANON;
		tcsetattr(STDIN, TCSANOW, &term);
		setbuf(stdin, 0);
		initialized = 1;
	}
	int bytes_waiting;
	ioctl(STDIN, FIONREAD, &bytes_waiting);
	return bytes_waiting;
}

int getch_toupper() {
	int button = toupper(__getch());
	if (button == '\033') {
		unsigned char new_button = __getch();
		new_button = __getch();
		if (new_button == 'A') { button = 'W'; }
		if (new_button == 'B') { button = 'S'; }
		if (new_button == 'D') { button = 'A'; }
		if (new_button == 'C') { button = 'D'; }
	}
	return button;
}

inline void print_console(const char* str, int length) {
	printf("%s", str); //std::cout << std::flush;
}
#endif // Linux

class ImageBmp {
public:
	~ImageBmp() {
		crop(0, 0);
	}

	ImageBmp(int width, int height) {
		crop(width, height);
	}

	inline void crop(int new_width, int new_height) {
		bmp.resize(new_height, std::vector<unsigned char>(new_width * 3));
		width = new_width;
		height = new_height;
	}

	inline void save_as(const std::string& filename) {
		stream.open(filename, std::ofstream::binary);
		int iter = 0;
		int extra_bytes = 4 - ((width * 3) % 4);
		extra_bytes = extra_bytes == 4 ? 0 : extra_bytes;
		unsigned int padded_size = ((width * 3) + extra_bytes) * height;
		unsigned int one_headers[6]{ padded_size + 54, 0, 54, 40, (unsigned int)width, (unsigned int)height };
		unsigned int two_headers[6]{ 0, padded_size, 0, 0, 0, 0 };
		std::vector<unsigned char> line_headers(54);
		line_headers[iter++] = 'B';
		line_headers[iter++] = 'M';
		for (int i = 0; i < 6; ++i) {
			line_headers[iter++] = (unsigned char)((one_headers[i] & 0x000000ff));
			line_headers[iter++] = (unsigned char)((one_headers[i] & 0x0000ff00) >> 8);
			line_headers[iter++] = (unsigned char)((one_headers[i] & 0x00ff0000) >> 16);
			line_headers[iter++] = (unsigned char)((one_headers[i] & (unsigned int)0xff000000) >> 24);
		}
		line_headers[iter++] = 1;
		line_headers[iter++] = 0;
		line_headers[iter++] = 24;
		line_headers[iter++] = 0;
		for (int i = 0; i < 6; ++i) {
			line_headers[iter++] = (unsigned char)((two_headers[i] & 0x000000ff));
			line_headers[iter++] = (unsigned char)((two_headers[i] & 0x0000ff00) >> 8);
			line_headers[iter++] = (unsigned char)((two_headers[i] & 0x00ff0000) >> 16);
			line_headers[iter++] = (unsigned char)((two_headers[i] & (unsigned int)0xff000000) >> 24);
		}
		stream.write(reinterpret_cast<const char*>(line_headers.data()), line_headers.size());
		if (!bmp[0].empty()) {
			std::streamsize stream_size = bmp[0].size();
			for (int y = height - 1; y >= 0; --y) {
				stream.write(reinterpret_cast<const char*>(bmp[y].data()), stream_size);
			}
		}
		stream.close();
	}

	inline void read_as(const std::string& filename) {
		std::ifstream stream(filename.c_str(), std::ifstream::binary);
		unsigned char line_headers[54];
		stream.read((char*)line_headers, 54);
		crop(*(int*)&line_headers[18], *(int*)&line_headers[22]);
		unsigned int offset = *(unsigned int*)&line_headers[10];
		stream.ignore(offset - 54);
		for (int y = height - 1; y >= 0; --y) {
			for (int x = 0; x < width; ++x) {
				int x_3 = x * 3;
				bmp[y][x_3] = stream.get();
				bmp[y][x_3 + 1] = stream.get();
				bmp[y][x_3 + 2] = stream.get();
			}
		}
		stream.close();
	}

	inline void set_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) noexcept {
		int x_3 = x * 3;
		bmp[y][x_3] = b;
		bmp[y][x_3 + 1] = g;
		bmp[y][x_3 + 2] = r;
	}

	inline int get_width() const noexcept {
		return width;
	}

	inline int get_height() const noexcept {
		return height;
	}

private:
	int width = 0;
	int height = 0;
	std::ofstream stream;
	std::vector<std::vector<unsigned char>> bmp;
};

inline void mandelbrot(double* arrays, int fractal_index, double z_x, double z_i, double c_x, double c_i, double factor, double power, double range, double inverse, double mult_power, int ismandel, int width, int height, int max_iteration, long max_sample) {
	double x_set = ismandel ? c_x : z_x;
	double i_set = ismandel ? c_i : z_i;
	double x_min = x_set - 1.5 * factor;
	double x_max = x_set + 1.5 * factor;
	double y_min = i_set - factor;
	double y_max = i_set + factor;
	double dx = (x_max - x_min) / (width - 1);
	double dy = (y_max - y_min) / (height - 1);
	double c_real, c_imag, z_real, z_imag, z_power, z_real_power = 0.0, z_imag_power = 0.0, z_power_sum = 0.0, z_atan2, power_sum_invert, z_imag_prev = 0.0, z_real_prev = 0.0;
	double z_cos_real, z_cos_imag, z_right_real, z_right_imag, PI = 3.14159265359, kx, ky, z_power_invert, smooth_value = 0.0, nu = 0.0;
	double z_temp_real = 0.0, z_temp_imag = 0.0, z_real2 = 0.0, z_imag2 = 0.0, p_real = -0.5, p_imag = 0.0, z_real_minus_tmp_real_power = 0.0, z_imag_minus_tmp_imag_power = 0.0;
	int n = 0, index = 0, escape_newton = 0;
	int is_power_two = power == 2.0 && mult_power == 1.0, is_not_power_one = !(power == 1.0 && mult_power == 1.0);
	double log_value = power > 2 ? 1.0 / fast_log(power > 0.0 ? power : -power) : 1.44269504088896, p_y = 0.0;
	double half_log_value = log_value * 0.5, real_power = power * 0.5 * mult_power, power_minus_one = power - 1.0, real_newton_power = (power - 1.0) * 0.5 * mult_power;
	if (fractal_index == 9) { range = pow(10, 21); }
	if (fractal_index == 6) { range = pow(range, 21); }
	if (fractal_index == 5) { range = pow(range, 12); }
	if (fractal_index == 4) { range = pow(range, 7); }
	if (fractal_index == 3) { range = pow(range, 2.618); }
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			c_real = ismandel ? x_min + x * dx : c_x;
			c_imag = ismandel ? y_min + y * dy : c_i;
			if (inverse != 0.0) {
				power_sum_invert = 1.0 / (c_real * c_real + c_imag * c_imag) * inverse;
				c_real = power_sum_invert * c_real + (1.0 - inverse) * c_real;
				c_imag = power_sum_invert * c_imag + (1.0 - inverse) * c_imag;
			}
			z_real = ismandel ? c_real : x_min + x * dx;
			z_imag = ismandel ? c_imag : y_min + y * dy;
			if (fractal_index == 5) { z_imag = -z_imag; }
			z_imag2 = 0.0, z_real2 = 0.0;
			n = 1;
			while (n < max_iteration) {
				z_real_power = z_real * z_real;
				z_imag_power = z_imag * z_imag;
				z_power_sum = z_real_power + z_imag_power;
				if (z_power_sum > range) {
					break;
				}
				if (fractal_index == 2) { // burning ship
					z_real = fabs(z_real);
					z_imag = fabs(z_imag);
				}
				if (is_power_two) { // z ^ 2
					z_imag_prev = z_imag * 2.0 * z_real;
					z_real_prev = z_real_power - z_imag_power;
				}
				else if (is_not_power_one) { // z ^ n
					z_power = fast_pow(z_power_sum, real_power);
					z_atan2 = atan2(z_imag, z_real);
					z_imag_prev = z_power * sin(power * z_atan2);
					z_real_prev = z_power * cos(power * z_atan2);
				}
				else {
					z_imag_prev = z_imag;
					z_real_prev = z_real;
				}
				if (fractal_index < 3) { // mandelbrot
					if (fractal_index == 1) { // tricorn
						z_imag_prev = -z_imag_prev;
					}
					z_real = z_real_prev;
					z_imag = z_imag_prev;
				}
				else {
					if (fractal_index == 3) { // collatz
						kx = exp(-z_imag_prev * PI);
						ky = exp(z_imag_prev * PI);
						z_cos_real = cos(z_real_prev * PI) * (0.5 * (ky + kx));
						z_cos_imag = sin(z_real_prev * PI) * (0.5 * (ky - kx));
						z_right_real = 5.0 * z_real_prev + 2.0;
						z_right_imag = 5.0 * z_imag_prev;
						z_real = (7.0 * z_real_prev + 2.0 - (z_cos_real * z_right_real + z_cos_imag * z_right_imag)) * 0.25;
						z_imag = (7.0 * z_imag_prev - (z_cos_real * z_right_imag - z_cos_imag * z_right_real)) * 0.25;
					}
					if (fractal_index == 4) { // collatz_mandelbrot
						kx = exp(z_real_prev);
						z_cos_real = cos(z_imag_prev) * kx;
						z_cos_imag = sin(z_imag_prev) * kx;
						z_right_real = 5.0 * z_real_prev + 2.0;
						z_right_imag = 5.0 * z_imag_prev;
						z_real = (7.0 * z_real_prev + 2.0 - (z_cos_real * z_right_real - z_cos_imag * z_right_imag)) * 0.25;
						z_imag = (7.0 * z_imag_prev - (z_cos_real * z_right_imag + z_cos_imag * z_right_real)) * 0.25;
					}
					if (fractal_index == 5) { // collatz_like_mandelbrot
						kx = exp(-z_imag_prev * PI);
						z_cos_real = cos(z_real_prev * PI) * kx;
						z_cos_imag = sin(z_real_prev * PI) * kx;
						z_right_real = 5.0 * z_real_prev + 2.0;
						z_right_imag = 5.0 * z_imag_prev;
						z_real = (7.0 * z_real_prev + 2.0 - (z_cos_real * z_right_real - z_cos_imag * z_right_imag)) * 0.25;
						z_imag = (7.0 * z_imag_prev - (z_cos_real * z_right_imag + z_cos_imag * z_right_real)) * 0.25;
					}
					if (fractal_index == 6) { // e ^ z
						kx = exp(z_real_prev);
						z_cos_real = cos(z_imag_prev) * kx;
						z_cos_imag = sin(z_imag_prev) * kx;
						z_real = z_cos_real;
						z_imag = z_cos_imag;
					}
					if (fractal_index == 7) { // phoenix
						z_temp_real = z_real;
						z_temp_imag = z_imag;
						z_real = z_real_prev + (z_real2 * p_real - z_imag2 * p_imag);
						z_imag = z_imag_prev + (z_real2 * p_imag + z_imag2 * p_real);
						z_real2 = z_temp_real;
						z_imag2 = z_temp_imag;
					}
					if (fractal_index == 8) { // feather
						z_real_power++;
						z_power_invert = 1.0 / (z_real_power * z_real_power + z_imag_power * z_imag_power);
						z_real = (z_real_prev * z_real_power + z_imag_prev * z_imag_power) * z_power_invert;
						z_imag = (z_imag_prev * z_real_power - z_real_prev * z_imag_power) * z_power_invert;
					}
					if (fractal_index == 9) { // newton
						z_power = power * fast_pow(z_power_sum, real_newton_power);
						if (is_power_two || !is_not_power_one) {
							z_atan2 = atan2(z_imag, z_real);
						}
						z_imag_power = z_power * sin(power_minus_one * z_atan2);
						z_real_power = z_power * cos(power_minus_one * z_atan2);
						z_temp_real = z_real;
						z_temp_imag = z_imag;
						z_real_prev--;
						z_power_invert = 1.0 / (z_real_power * z_real_power + z_imag_power * z_imag_power);
						z_real = z_real - (z_real_prev * z_real_power + z_imag_prev * z_imag_power) * z_power_invert;
						z_imag = z_imag - (z_imag_prev * z_real_power - z_real_prev * z_imag_power) * z_power_invert;
						z_real_minus_tmp_real_power = z_real + c_real - z_temp_real;
						z_imag_minus_tmp_imag_power = z_imag + c_imag - z_temp_imag;
						if (z_real_minus_tmp_real_power * z_real_minus_tmp_real_power + z_imag_minus_tmp_imag_power * z_imag_minus_tmp_imag_power < 0.002) {
							break;
						}
					}
				}
				z_real += c_real;
				z_imag += c_imag; // z_temp_imag
				n++;
			}
			smooth_value = 0.0;
			if (n != 1 && n != max_iteration && fractal_index != 9) {
				nu = fast_log(fast_log(z_power_sum) * half_log_value) * log_value;
				if (nu < n) { smooth_value = (double)n - nu; }
			}
			else if (fractal_index == 9) {
				smooth_value = n == max_iteration ? 0 : n;
			}
			arrays[index++] = smooth_value;
		}
	}
}

inline void mandelbulb(double* arrays, int fractal_index, double z_x, double z_y, double cx, double cy, double factor, double power, double range, double inverse, double mult_power, int ismandel, int width, int height, int max_iteration, long max_sample, double min_dist, double rot_x, double rot_z, int raymarch_iterations, double phi_shift, int formula, double cz, double cw, double zz, double zw) {
	double x_set = ismandel ? cx : z_x;
	double i_set = ismandel ? cy : z_y;
	x_set = i_set = 0;
	double factor_one = 1.0;
	double x_min = x_set - 1.5 * factor_one;
	double x_max = x_set + 1.5 * factor_one;
	double y_min = i_set - factor_one;
	double y_max = i_set + factor_one;
	double dx = (x_max - x_min) / (width - 1);
	double dy = (y_max - y_min) / (height - 1);
	double cam_ro_x = (ismandel ? cx : z_x);
	double cam_ro_y = (ismandel ? cy : z_y);
	double cam_ro_z = -factor;
	double real_power = power * mult_power, surf_dist = 0.0, z_x_tmp, z_y_tmp, z_z_tmp, z_x_2, z_y_2, z_z_2, p_x = -0.5, p_y = 0.0, p_z = 0.0;
	double cos_rot_x = cos(rot_x), sin_rot_x = sin(rot_x), cos_rot_z = cos(rot_z), sin_rot_z = sin(rot_z), inv_raymarch_iterations = 1.0 / (double)raymarch_iterations;
	int is_power_two = power == 2.0 && phi_shift == 0.0 && mult_power == 1.0, max_iterations = max_iteration / 5, index = 0;
	range = range * 4.0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double uv_x = x_min + x * dx;
			double uv_y = y_min + y * dy;
			double rd_z = 1.0 / sqrt(uv_x * uv_x + uv_y * uv_y + 1.0);
			double rd_y = uv_y * rd_z;
			double rd_x = uv_x * rd_z;
			double col = 0.0, dO = 0.0;
			int ot = raymarch_iterations - 1;
			for (int i = 0; i < raymarch_iterations; i++) {
				double cam_z_x = cam_ro_x + rd_x * dO;
				double cam_z_y = cam_ro_y + rd_y * dO;
				double cam_z_z = cam_ro_z + rd_z * dO;
				vec3 z = vec3{ cam_z_x, cam_z_y, cam_z_z };
				double c = cos_rot_x;
				double s = sin_rot_x;
				double tmp_z = c * z.z - s * z.y;
				z.y = c * z.y + s * z.z;
				z.z = tmp_z;
				c = cos_rot_z;
				s = sin_rot_z;
				tmp_z = c * z.z - s * z.x;
				z.x = c * z.x + s * z.z;
				z.z = tmp_z;
				double dr = 1.0, a = 0.0, r = 0.0, sin_theta = 0.0, zr, phi, theta, cos_phi, dist = 0.0;
				double z_x_prev = 0.0, z_y_prev = 0.0, z_z_prev = 0.0;
				int orbit = 0;
				double c_x = ismandel ? z.x : cx;
				double c_y = ismandel ? z.y : cy;
				double c_z = ismandel ? z.z : 0.0;
				if (inverse != 0.0) {
					double c_pow_sum = inverse * 1.0 / (c_x * c_x + c_y * c_y + c_z * c_z);
					c_x = c_x * c_pow_sum + ((1.0 - inverse) * c_x);
					c_y = c_y * c_pow_sum + ((1.0 - inverse) * c_y);
					c_z = c_z * c_pow_sum + ((1.0 - inverse) * c_z);
				}
				double z_x = ismandel ? c_x : z.x;
				double z_y = ismandel ? c_y : z.y;
				double z_z = ismandel ? c_z : z.z;
				for (int j = 1; j < max_iterations; j++) {
					double zx2 = z_x * z_x;
					double zy2 = z_y * z_y;
					double zz2 = z_z * z_z;
					double zxy2 = zx2 + zy2;
					r = sqrt(zxy2 + zz2);
					if (r > range) {
						orbit = j;
						break;
					}
					if (fractal_index == 2) {
						z_x = fabs(z_x);
						z_y = fabs(z_y);
						z_z = -fabs(z_z);
					}
					if (formula == 0) {
						if (is_power_two) {
							dr = r * dr * 2.0 + 1.0;
							a = 1.0 - zz2 / zxy2;
							z_y_prev = 2.0 * z_y * z_x * a;
							z_z_prev = 2.0 * z_z * sqrt(zxy2);
							z_x_prev = (zx2 - zy2) * a;
						}
						else {
							dr = fast_pow(r, power - 1.0) * dr * power + 1.0;
							phi = asin(z_z / r) + phi_shift;
							theta = atan2(z_y, z_x);
							zr = fast_pow(r, real_power);
							cos_phi = cos(power * phi);
							z_x_prev = zr * cos(power * theta) * cos_phi;
							z_y_prev = zr * sin(power * theta) * cos_phi;
							z_z_prev = zr * sin(power * phi);
						}
					}
					else if (formula == 1) {
						zr = pow(r, power);
						theta = acos(z_z / r);
						phi = atan2(z_y, z_x);
						sin_theta = sin(power * theta);
						z_x_prev = sin_theta * cos(power * phi) * zr;
						z_y_prev = sin_theta * sin(power * phi) * zr;
						z_z_prev = cos(power * theta) * zr;
						dr = pow(r, power - 1.0) * dr * power + 1.0;
					}
					if (fractal_index < 3) {
						if (fractal_index == 1) { // tricorn
							z_y_prev = -z_y_prev;
							z_z_prev = z_z_prev; // minus
						}
						z_x = z_x_prev;
						z_y = z_y_prev;
						z_z = z_z_prev;
					}
					z_x = z_x + c_x;
					z_y = z_y + c_y;
					z_z = z_z + c_z;
				}
				dist = 0.5 * log(r) * r / dr;
				dO += dist;
				if (surf_dist == 0.0) { surf_dist = dist * min_dist; }
				if (dist < surf_dist) { ot = i; col = (double)orbit; break; }
				if (dO > 100.0) { ot = 0; col = 0; break; }
			}
			double together = 0.0;
			if (ot != 0) {
				double ot_div = (double)ot * inv_raymarch_iterations;
				together = col + (ot_div >= 1.0 ? ot_div - 0.001 : ot_div);
			}
			arrays[index++] = together; // 13.909 | 13 - it's col value, 909 - it's ot_div value
		}
	}
}

inline void mandel4d(double* arrays, int fractal_index, double z_x, double z_y, double cx, double cy, double factor, double power, double range, double inverse, double mult_power, int ismandel, int width, int height, int max_iteration, long max_sample, double min_dist, double rot_x, double rot_z, int raymarch_iterations, double phi_shift, int formula, double cz, double cw, double zz, double zw) {
	double x_set = ismandel ? cx : z_x;
	double i_set = ismandel ? cy : z_y;
	x_set = i_set = 0;
	double factor_one = 1.0;
	double x_min = x_set - 1.5 * factor_one;
	double x_max = x_set + 1.5 * factor_one;
	double y_min = i_set - factor_one;
	double y_max = i_set + factor_one;
	double dx = (x_max - x_min) / (width - 1);
	double dy = (y_max - y_min) / (height - 1);
	double cam_ro_x = (ismandel ? cx : z_x);
	double cam_ro_y = (ismandel ? cy : z_y);
	double cam_ro_z = -factor;
	double real_power = power * mult_power, surf_dist = 0.0, z_x_tmp, z_y_tmp, z_z_tmp, z_x_2, z_y_2, z_z_2, p_x = -0.5, p_y = 0.0, p_z = 0.0;
	double cos_rot_x = cos(rot_x), sin_rot_x = sin(rot_x), cos_rot_z = cos(rot_z), sin_rot_z = sin(rot_z), inv_raymarch_iterations = 1.0 / (double)raymarch_iterations;
	int is_power_two = power == 2.0 && phi_shift == 0.0 && mult_power == 1.0, max_iterations = max_iteration, index = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double uv_x = x_min + x * dx;
			double uv_y = y_min + y * dy;
			double rd_z = 1.0 / sqrt(uv_x * uv_x + uv_y * uv_y + 1.0);
			double rd_y = uv_y * rd_z;
			double rd_x = uv_x * rd_z;
			double col = 0.0, dO = 0.0;
			int ot = raymarch_iterations - 1;
			for (int i = 0; i < raymarch_iterations; i++) {
				double cam_z_x = cam_ro_x + rd_x * dO;
				double cam_z_y = cam_ro_y + rd_y * dO;
				double cam_z_z = cam_ro_z + rd_z * dO;
				vec3 z = vec3{ cam_z_x, cam_z_y, cam_z_z };
				double c = cos_rot_x;
				double s = sin_rot_x;
				double tmp_z = c * z.z - s * z.y;
				z.y = c * z.y + s * z.z;
				z.z = tmp_z;
				c = cos_rot_z;
				s = sin_rot_z;
				tmp_z = c * z.z - s * z.x;
				z.x = c * z.x + s * z.z;
				z.z = tmp_z;
				double dr = 1.0, a = 0.0, r = 0.0, sin_theta = 0.0, zr, phi, theta, cos_phi, dist = 0.0;
				double z_x_prev = 0.0, z_y_prev = 0.0, z_z_prev = 0.0, z_w_prev = 0.0;
				int orbit = 0;
				double c_x = ismandel ? z.x : cx;
				double c_y = ismandel ? z.y : cy;
				double c_z = ismandel ? z.z : cz;
				double c_w = ismandel ? zw : cw;
				if (inverse != 0.0) {
					double c_pow_sum = inverse * 1.0 / (c_x * c_x + c_y * c_y + c_z * c_z + c_w * c_w);
					c_x = c_x * c_pow_sum + ((1.0 - inverse) * c_x);
					c_y = c_y * c_pow_sum + ((1.0 - inverse) * c_y);
					c_z = c_z * c_pow_sum + ((1.0 - inverse) * c_z);
					c_w = c_w * c_pow_sum + ((1.0 - inverse) * c_w);
				}
				double z_x = ismandel ? c_x : z.x;
				double z_y = ismandel ? c_y : z.y;
				double z_z = ismandel ? c_z : z.z;
				double z_w = ismandel ? cw : zw;
				for (int j = 1; j < max_iterations / 5; j++) {
					double zx2 = z_x * z_x;
					double zy2 = z_y * z_y;
					double zz2 = z_z * z_z;
					double zw2 = z_w * z_w;
					r = sqrt(zx2 + zy2 + zz2 + zw2);
					if (r > range) {
						orbit = j;
						break;
					}
					if (fractal_index == 2) {
						z_x = fabs(z_x);
						z_y = fabs(z_y);
						z_z = fabs(z_z);
						z_w = fabs(z_w);
					}
					if (is_power_two) {
						dr = fast_pow(r, power - 1.0) * dr * power + 1.0;
						if (formula == 0) { // quaternion
							z_x_prev = zx2 - zy2 - zz2 - zw2;
							z_y_prev = 2.0 * z_x * z_y;
							z_z_prev = 2.0 * z_x * z_z;
							z_w_prev = 2.0 * z_x * z_w;
						}
						if (formula == 1) { // bicomplex
							z_x_prev = zx2 - zy2 - zz2 + zw2;
							z_y_prev = 2.0 * (z_x * z_y - z_z * z_w);
							z_z_prev = 2.0 * (z_x * z_z - z_y * z_w);
							z_w_prev = 2.0 * (z_x * z_w + z_y * z_z);
						}
					}
					if (fractal_index < 3) {
						if (fractal_index == 1) {
							z_y_prev = -z_y_prev;
							z_z_prev = -z_z_prev;
							z_w_prev = -z_w_prev;
						}
						z_x = z_x_prev;
						z_y = z_y_prev;
						z_z = z_z_prev;
						z_w = z_w_prev;
					}
					z_x = z_x + c_x;
					z_y = z_y + c_y;
					z_z = z_z + c_z;
					z_w = z_w + c_w;
				}
				dist = 0.5 * log(r) * r / dr;
				dO += dist;
				if (surf_dist == 0.0) { surf_dist = dist * min_dist; }
				if (dist < surf_dist) { ot = i; col = (double)orbit; break; }
				if (dO > 100.0) { ot = 0; col = 0; break; }
			}
			double together = 0.0;
			if (ot != 0) {
				double ot_div = (double)ot * inv_raymarch_iterations;
				together = col + (ot_div >= 1.0 ? ot_div - 0.001 : ot_div);
			}
			arrays[index++] = together; // 13.909 | 13 - it's col value, 909 - it's ot_div value
		}
	}
}

unsigned int current = 1;
inline void buddhabrot(double* arrays, int fractal_index, double z_x, double z_i, double c_x, double c_i, double factor, double power, double range, double inverse, double mult_power, int ismandel, int width, int height, int iteration, long max_sample) {
	double x_set = ismandel ? c_x : z_x;
	double i_set = ismandel ? c_i : z_i;
	double x_min = x_set - 1.5 * factor;
	double x_max = x_set + 1.5 * factor;
	double y_min = i_set - factor;
	double y_max = i_set + factor;
	double x_size = x_max - x_min;
	double y_size = y_max - y_min;
	double nx_factor = 1.0 / x_size * width;
	double ny_factor = 1.0 / y_size * height;
	double c_real, c_imag, z_real, z_imag, z_power_sum, z_real_power, z_imag_power, z_atan2, power_sum_invert, z_imag_prev, z_real_prev;
	double z_cos_real, z_cos_imag, z_right_real, z_right_imag, PI = 3.14159265359, kx, ky, z_power_invert, z_power;
	double z_temp_real = 0.0, z_temp_imag = 0.0, z_old_real = 0.0, z_old_imag = 0.0, z_real2 = 0.0, z_imag2 = 0.0, p_real = -0.5, p_imag = 0.0, z_real_minus_tmp_real_power = 0.0, z_imag_minus_tmp_imag_power = 0.0;
	double rand_real, rand_imag, mult_distrub = 1.0 / (double)0xFFFFFFFF;
	double real_power = power * 0.5 * mult_power, power_minus_one = power - 1.0, real_newton_power = (power - 1.0) * 0.5 * mult_power;
	int is_power_two = power == 2.0 && mult_power == 1.0, is_power_one = power == 1.0 && mult_power == 1.0;
	int nx, ny, idx, area = width * height, escape_newton = 0, is_escape = 0;
	if (fractal_index == 6) { range = pow(range, 21.0); }
	for (long i = 0; i < max_sample; i++) {
		z_real2 = z_imag2 = 0.0;
		current = (166421U * current + 1054222352U) % 0xFFFFFFFF;
		rand_real = current * mult_distrub;
		current = (1664525U * current + 1013904223U) % 0xFFFFFFFF;
		rand_imag = current * mult_distrub;
		c_real = ismandel ? (rand_real * x_size + x_min) : c_x;
		c_imag = ismandel ? (rand_imag * y_size + y_min) : c_i;
		if (inverse != 0.0) {
			power_sum_invert = 1.0 / (c_real * c_real + c_imag * c_imag) * inverse;
			c_real = c_real * power_sum_invert + ((1.0 - inverse) * c_real);
			c_imag = c_imag * power_sum_invert + ((1.0 - inverse) * c_imag);
		}
		z_real = ismandel ? c_real : (rand_real * x_size + x_min);
		z_imag = ismandel ? c_imag : (rand_imag * y_size + y_min);
		is_escape = 0; escape_newton = 0;
		z_old_real = z_real;
		z_old_imag = z_imag;
		for (int k = 0; k < iteration; k++) {
			z_real_power = z_real * z_real;
			z_imag_power = z_imag * z_imag;
			z_power_sum = z_real_power + z_imag_power;
			if (fractal_index == 2) { // burning ship
				z_real = fabs(z_real);
				z_imag = fabs(z_imag);
			}
			if (is_power_two) { // z ^ 2
				z_imag_prev = 2.0 * z_imag * z_real;
				z_real_prev = z_real_power - z_imag_power;
			}
			else if (is_power_one) {
				z_imag_prev = z_imag;
				z_real_prev = z_real;
			}
			else { // z ^ n
				z_power = fast_pow(z_power_sum, real_power);
				z_atan2 = atan2(z_imag, z_real);
				z_imag_prev = z_power * sin(power * z_atan2);
				z_real_prev = z_power * cos(power * z_atan2);
			}
			if (fractal_index < 3) { // mandelbrot
				if (fractal_index == 1) { // tricorn
					z_imag_prev = -z_imag_prev;
				}
				z_real = z_real_prev + c_real;
				z_imag = z_imag_prev + c_imag;
			}
			if (fractal_index == 3) { // collatz
				kx = exp(-z_imag_prev * PI);
				ky = exp(z_imag_prev * PI);
				z_cos_real = cos(z_real_prev * PI) * (0.5 * (ky + kx));
				z_cos_imag = sin(z_real_prev * PI) * (0.5 * (ky - kx));
				z_right_real = 5.0 * z_real_prev + 2.0;
				z_right_imag = 5.0 * z_imag_prev;
				z_real = (7.0 * z_real_prev + 2.0 - (z_cos_real * z_right_real + z_cos_imag * z_right_imag)) * 0.25 + c_real;
				z_imag = (7.0 * z_imag_prev - (z_cos_real * z_right_imag - z_cos_imag * z_right_real)) * 0.25 + c_imag;
			}
			if (fractal_index == 4) { // collatz_mandelbrot
				kx = exp(z_real_prev);
				z_cos_real = cos(z_imag_prev) * kx;
				z_cos_imag = sin(z_imag_prev) * kx;
				z_right_real = 5.0 * z_real_prev + 2.0;
				z_right_imag = 5.0 * z_imag_prev;
				z_real = (7.0 * z_real_prev + 2.0 - (z_cos_real * z_right_real - z_cos_imag * z_right_imag)) * 0.25 + c_real;
				z_imag = (7.0 * z_imag_prev - (z_cos_real * z_right_imag + z_cos_imag * z_right_real)) * 0.25 + c_imag;
			}
			if (fractal_index == 5) { // collatz_like_mandelbrot
				kx = exp(-z_imag_prev * PI);
				z_cos_real = cos(z_real_prev * PI) * kx;
				z_cos_imag = sin(z_real_prev * PI) * kx;
				z_right_real = 5.0 * z_real_prev + 2.0;
				z_right_imag = 5.0 * z_imag_prev;
				z_real = (7.0 * z_real_prev + 2.0 - (z_cos_real * z_right_real - z_cos_imag * z_right_imag)) * 0.25 + c_real;
				z_imag = (7.0 * z_imag_prev - (z_cos_real * z_right_imag + z_cos_imag * z_right_real)) * 0.25 + c_imag;
			}
			if (fractal_index == 6) { // e ^ z
				kx = exp(z_real_prev);
				z_cos_real = cos(z_imag_prev) * kx;
				z_cos_imag = sin(z_imag_prev) * kx;
				z_real = z_cos_real + c_real;
				z_imag = z_cos_imag + c_imag;
			}
			if (fractal_index == 7) { // phoenix
				z_temp_real = z_real;
				z_temp_imag = z_imag;
				z_real = z_real_prev + (z_real2 * p_real - z_imag2 * p_imag) + c_real;
				z_imag = z_imag_prev + (z_real2 * p_imag + z_imag2 * p_real) + c_imag;
				z_real2 = z_temp_real;
				z_imag2 = z_temp_imag;
			}
			if (fractal_index == 8) { // feather
				z_real_power++;
				z_power_invert = 1.0 / (z_real_power * z_real_power + z_imag_power * z_imag_power);
				z_real = (z_real_prev * z_real_power + z_imag_prev * z_imag_power) * z_power_invert + c_real;
				z_imag = (z_imag_prev * z_real_power - z_real_prev * z_imag_power) * z_power_invert + c_imag;
			}
			if (fractal_index == 9) { // newton
				z_real_prev--;
				z_power = power * fast_pow(z_power_sum, real_newton_power);
				if (is_power_two || is_power_one) {
					z_atan2 = atan2(z_imag, z_real);
				}
				z_imag_power = z_power * sin(power_minus_one * z_atan2);
				z_real_power = z_power * cos(power_minus_one * z_atan2);
				z_temp_real = z_real;
				z_temp_imag = z_imag;
				z_power_invert = 1.0 / (z_real_power * z_real_power + z_imag_power * z_imag_power);
				z_real = (z_real - (z_real_power * z_real_prev + z_imag_power * z_imag_prev) * z_power_invert + c_real);
				z_imag = (z_imag - (z_real_power * z_imag_prev - z_imag_power * z_real_prev) * z_power_invert + c_imag);
				z_real_minus_tmp_real_power = z_real - z_temp_real;
				z_imag_minus_tmp_imag_power = z_imag - z_temp_imag;
				escape_newton = z_real_minus_tmp_real_power * z_real_minus_tmp_real_power + z_imag_minus_tmp_imag_power * z_imag_minus_tmp_imag_power < 0.00275 && fractal_index == 9;
			}
			if (is_escape) {
				if (z_power_sum > range && fractal_index != 9 || escape_newton) {
					break;
				}
				if (z_real >= x_max || z_real < x_min || z_imag >= y_max || z_imag < y_min) {
					continue;
				}
				nx = (int)((z_real - x_min) * nx_factor);
				ny = (int)((z_imag - y_min) * ny_factor);
				idx = nx + ny * width;
				if (idx < area) {
					arrays[idx]++;
				}
			}
			if (z_power_sum > range && fractal_index != 9 || escape_newton) {
				is_escape = 1;
				z_real = z_old_real;
				z_imag = z_old_imag;
				z_real2 = z_imag2 = 0.0;
			}
		}
	}
}

inline void lyapunov(double* arrays, double z_x, double z_i, double c_x, double c_i, double factor, double power, double range, double inverse, double mult_power, int ismandel, int width, int height, int max_iteration, long max_sample) {
	double x_set = ismandel ? c_x : z_x;
	double i_set = ismandel ? c_i : z_i;
	double x_min = x_set - 1.5 * factor;
	double x_max = x_set + 1.5 * factor;
	double y_min = i_set - factor;
	double y_max = i_set + factor;
	double dx = (x_max - x_min) / (width - 1);
	double dy = (y_max - y_min) / (height - 1);
	double c_real, c_imag, z_real, z_imag, power_sum_invert;
	double lambda, max_lambda = 120000.0, xn = 0.5, r = 0.0, ri = 0.0, iteration_less = 1.0 / (double)max_iteration;
	int coord_use[] = { 1, 0 };
	int coord_count = (int)sizeof(coord_use) / sizeof(coord_use[0]), value = 0;
	int index = 0, n0 = clamp(min(max_iteration / 5, max_iteration - max_iteration / 5), 1, max_iteration);
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			c_real = ismandel ? x_min + x * dx : c_x;
			c_imag = ismandel ? y_min + y * dy : c_i;
			if (inverse != 0.0) {
				power_sum_invert = 1.0 / (c_real * c_real + c_imag * c_imag) * inverse;
				c_real = power_sum_invert * c_real + (1.0 - inverse) * c_real;
				c_imag = power_sum_invert * c_imag + (1.0 - inverse) * c_imag;
			}
			z_real = ismandel ? c_real : x_min + x * dx;
			z_imag = -(ismandel ? c_imag : y_min + y * dy); // - minus
			lambda = ri = 0.0;
			xn = 0.5;
			for (int n = 1; n < max_iteration; n++) {
				switch (coord_use[n % coord_count]) {
				case 0: r = z_real; ri = c_imag; break;
				case 1: r = z_imag; ri = c_real; break;
				}
				xn = xn * r * (1.0 - xn);
				if (n >= n0) {
					lambda += log(fabs(r * (1.0 - 2.0 * xn))) * iteration_less;
				}
				if (fabs(lambda) > max_lambda) {
					break;
				}
			}
			if (lambda > 0.0) {
				value = (int)(exp(-3.3 * lambda) * 25.5);
			}
			else {
				lambda = exp(lambda);
				value = (int)(lambda <= 0.98 ? lambda * 15.2653061226 : (1.76 + -37.0 * (0.98 - lambda)) * 8.5);
			}
			arrays[index++] = value;
		}
	}
}

inline void get_fractal(double* arrays, int fractal_index, int algorithm, double z_x, double z_i, double c_x, double c_i, double factor, double power, double range, double inverse, double mult_power, int ismandel, int width, int height, int iteration, long max_sample, double min_dist, double rot_x, double rot_y, double raymarch_iterations, double phi_shift, int formula, double cz, double cw, double zz, double zw) {
	if (power == 0.0 && fractal_index != 10 && inverse == 0) { for (int i = 0; i < height * width; i++) { arrays[i] = 0.0; } return; }
	if (algorithm == 0) { // escape algoritm (classic)
		if (fractal_index != 10) { mandelbrot(arrays, fractal_index, z_x, z_i, c_x, c_i, factor, power, range, inverse, mult_power, ismandel, width, height, iteration, max_sample); }
		if (fractal_index == 10) { lyapunov(arrays, z_x, z_i, c_x, c_i, factor, power, range, inverse, mult_power, ismandel, width, height, iteration, max_sample); }
	}
	if (algorithm == 1) { // buddha algorith (all bounce before escape)
		for (int i = 0; i < height * width; i++) { arrays[i] = 0.0; }
		if (fractal_index != 10) { buddhabrot(arrays, fractal_index, z_x, z_i, c_x, c_i, factor, power, range, inverse, mult_power, ismandel, width, height, iteration, max_sample); }
		if (fractal_index == 10) { /* don't know */ }
	}
	if (algorithm == 2) { // 3d (classic)
		mandelbulb(arrays, fractal_index, z_x, z_i, c_x, c_i, factor, power, range, inverse, mult_power, ismandel, width, height, iteration, max_sample, min_dist, rot_x, rot_y, raymarch_iterations, phi_shift, formula, cz, cw, zz, zw);
	}
	if (algorithm == 3) { // 4d (quaternion, bicomplex)
		mandel4d(arrays, fractal_index, z_x, z_i, c_x, c_i, factor, power, range, inverse, mult_power, ismandel, width, height, iteration, max_sample, min_dist, rot_x, rot_y, raymarch_iterations, phi_shift, formula, cz, cw, zz, zw);
	}
}

inline vec3 hue(double inp_x, double inp_y, double inp_z, double H) {
	double U = cos(H * 0.0174532925199);
	double W = sin(H * 0.0174532925199);
	vec3 ret;
	ret.x = (0.299 + 0.701 * U + 0.168 * W) * inp_x + (0.587 - 0.587 * U + 0.330 * W) * inp_y + (0.114 - 0.114 * U - 0.497 * W) * inp_z;
	ret.y = (0.299 - 0.299 * U - 0.328 * W) * inp_x + (0.587 + 0.413 * U + 0.035 * W) * inp_y + (0.114 - 0.114 * U + 0.292 * W) * inp_z;
	ret.z = (0.299 - 0.3 * U + 1.25 * W) * inp_x + (0.587 - 0.588 * U - 1.05 * W) * inp_y + (0.114 + 0.886 * U - 0.203 * W) * inp_z;
	return ret;
}

int main() {
	double z_x = 0.0, z_y = 0.0, z_z = 0.0, z_w = 0.0, c_x = 0.0, c_y = 0.0, c_z = 0.0, c_w = 0.0;
	double power = 2.0, range = 8.0, inverse = 0, mult_power = 1.0, scale_mult_power = 1.0, scale_inverse = 1.0, scale_power = 1.0;
	double factor_mandelbrot = 1.5, factor_julia = 1.5;
	double mult_sample = 8.0, mult_light = 0.0, light = 1.0 * 1.33 * 1.33, val = 0.0;
	double rot_x = 0.0, rot_z = 0.0, move_scale = 0.05, min_dist = 0.001, phi_shift = 0.0, phi_speed_mult = 1.0; // 3d
	int raymarch_iterations = 120, fomula_3d = 0, is_mandel4d = 0; // 3d/4d
	int width = 0, height = 0, old_width = 0, old_height = 0, render_width = 1920, render_height = 1080;
	int is_mandel = 1, is_buddha = 0, is_render = 0, is_ansi = use_ansi, is_delete_button = 0, is_mandelbulb = 0;
	int frame = 0, iteration = 100, type_fractal = 0, fractal_count = 11, fractal_index = 0, array_index = 0;
	int r = -2, g = -2, b = -2, tmp_r = -1, tmp_g = -1, tmp_b = -1;
	char key = '\b';
	double* prewiew = (double*)malloc(sizeof(double));
	double* prewiew_mandelbrot = (double*)malloc(sizeof(double));
	double* prewiew_julia = (double*)malloc(sizeof(double));
	double* render = (double*)malloc(render_width * render_height * sizeof(double));
	std::vector<char> repeats_button;
	std::string row, row1;
	ColorBit last_color = Black;
	ColorBit temp_color = Black;
	ColorBit style[] = { Blue, DarkCyan, Cyan, DarkGreen, Green, Yellow, DarkGray, Red, DarkRed, Black };
	char gradient[] = "..54SUEaPXAhHmO6bdqp9RGMDW#N8B%0$Q&";
	int gradient_size = sizeof(gradient) / sizeof(gradient[0]);
	int style_size = sizeof(style) / sizeof(style[0]);
	ImageBmp image(render_width, render_height);
	std::vector<std::vector<char>> push_button, sequence_button;
	while (1) {
		Size resolution = get_size_screen();
		width = resolution.width;
		height = resolution.height;
		gotoxy(0, 0);
		if (repeats_button.size() > 0 || key != '\b' || (width != old_width || height != old_height)) {
			if (width != old_width || height != old_height) {
				free(prewiew_julia);
				free(prewiew_mandelbrot);
				free(prewiew);
				prewiew = (double*)malloc(width * height * sizeof(double));
				prewiew_mandelbrot = (double*)malloc((width / 2) * height * sizeof(double));
				prewiew_julia = (double*)malloc((width / 2) * height * sizeof(double));
			}
			old_width = width;
			old_height = height;
			frame += 1;
			if (type_fractal == 0 || type_fractal == 1) {
				get_fractal(prewiew, fractal_index, is_buddha + is_mandelbulb * 2 + is_mandel4d * 3, z_x, z_y, c_x, c_y, is_mandel ? factor_mandelbrot : factor_julia, power, range, inverse, mult_power, is_mandel, width, height, is_buddha ? iteration / 1.25 : iteration, (width >= 200 ? mult_sample * 0.5 : mult_sample * 0.8) * width * height, min_dist, rot_x, rot_z, raymarch_iterations, phi_shift, fomula_3d, c_z, c_w, z_z, z_w);
			}
			if (type_fractal == 2 || type_fractal == 3) {
				get_fractal(prewiew_mandelbrot, fractal_index, is_buddha + is_mandelbulb * 2 + is_mandel4d * 3, z_x, z_y, c_x, c_y, factor_julia, power, range, inverse, mult_power, 0, width / 2, height, is_buddha ? iteration / 1.25 : iteration, (width >= 200 ? mult_sample * 0.5 : mult_sample * 0.8) * width * height, min_dist, rot_x, rot_z, raymarch_iterations, phi_shift, fomula_3d, c_z, c_w, z_z, z_w);
				get_fractal(prewiew_julia, fractal_index, is_buddha + is_mandelbulb * 2 + is_mandel4d * 3, z_x, z_y, c_x, c_y, factor_mandelbrot, power, range, inverse, mult_power, 1, width / 2, height, is_buddha ? iteration / 1.25 : iteration, (width >= 200 ? mult_sample * 0.5 : mult_sample * 0.8) * width * height, min_dist, rot_x, rot_z, raymarch_iterations, phi_shift, fomula_3d, c_z, c_w, z_z, z_w);
				int index = 0;
				for (int y = 0; y < height; y++) {
					if (type_fractal == 2) {
						for (int x = 0; x < width / 2; x++) {
							prewiew[index++] = prewiew_julia[x + y * (width / 2)];
						}
						for (int x = 0; x < width / 2; x++) {
							prewiew[index++] = prewiew_mandelbrot[x + y * (width / 2)];
						}
						if (width % 2 == 1) { prewiew[index++] = 0.0; }
					}
					if (type_fractal == 3) {
						for (int x = 0; x < width / 2; x++) {
							prewiew[index++] = prewiew_julia[x + y * (width / 2)];
							prewiew[index++] = prewiew_mandelbrot[x + y * (width / 2)];
						}
						if (width % 2 == 1) { prewiew[index++] = 0.0; }
					}
				}
			}
			array_index = 0;
			if (is_ansi) {
				mult_light = light * 1.3333333; // 1.3333333 | 1.15
				for (int y = 0; y < height - 1; y++) {
					for (int x = 0; x < width; x++) {
						val = prewiew[++array_index];
						if (is_mandelbulb || is_mandel4d) {
							double col = val <= 0 ? 0 : (int)val;
							double ot_div = val - col;
							r = g = b = 0;
							if (ot_div != 0.0) {
								double bright = (1.0 - 1.5 * pow(ot_div, 1.1)) * mult_light * 108.118042163;
								vec3 color = hue(0.85, 0.45, 0.15, 1.0 + clamp(col, 3.0, 20.0) * 80.0);
								r = clamp(color.x * bright, 0, 255);
								g = clamp(color.y * bright, 0, 255);
								b = clamp(color.z * bright, 0, 255);
							}
						}
						else {
							val = floor(val * 5.0) * 0.2 * mult_light;
							r = bounce(val, 0, 255, 5);
							g = bounce(val * 2.0, 0, 255, 5);
							b = bounce(val * 3.0, 0, 255, 5);
						}
						row.push_back('.');
						if (tmp_r != r || tmp_g != g || tmp_b != b) {
							row.append("\x1b[38;2;");
							if (r >= 100) { row.push_back(r / 100 + '0'); }
							if (r >= 10) { row.push_back(r / 10 % 10 + '0'); }
							row.push_back(r % 10 + '0');
							row.push_back(';');
							if (g >= 100) { row.push_back(g / 100 + '0'); }
							if (g >= 10) { row.push_back(g / 10 % 10 + '0'); }
							row.push_back(g % 10 + '0');
							row.push_back(';');
							if (b >= 100) { row.push_back(b / 100 + '0'); }
							if (b >= 10) { row.push_back(b / 10 % 10 + '0'); }
							row.push_back(b % 10 + '0');
							row.append(";7m");
							tmp_r = r; tmp_g = g; tmp_b = b;
						}
					}
				}
				print_console(row.c_str(), row.size());
				row.clear();
			}
			if (!is_ansi) {
				color(temp_color, Black);
				double iteration_less = 1.0 / (double)iteration * (style_size - 1);
				for (int y = 0; y < height - 1; y++) {
					for (int x = 0; x < width; x++) {
						val = prewiew[array_index++];
						if (is_mandelbulb || is_mandel4d) {
							double col = val <= 0 ? 0 : (int)val;
							double ot_div = val - col;
							if (ot_div != 0.0) {
								double bright = (1.0 - 1.5 * pow(ot_div, 1.1)) * 255.0;
								vec3 color = hue(0.15, 0.45, 0.85, 1.0 + clamp(col, 3, 20) * 80.0);
								val = (clamp(color.x * bright, 0, 255) + clamp(color.y * bright, 0, 255) + clamp(color.z * bright, 0, 255)) / 3.0;
							}
						}
						double value = bounce(val, 0, 255, 5) * 3.0;
						temp_color = style[(int)(value * iteration_less) % style_size];
						if (last_color != temp_color) {
							print_console(row1.c_str(), row1.size());
							row1.clear();
							last_color = temp_color;
							color(temp_color, Black);
						}
						row1.push_back(gradient[(int)value % gradient_size]);
					}
					row1.push_back('\n');
					print_console(row1.c_str(), row1.size());
					row1.clear();
				}
			}
			gotoxy(0, 0);
			if (is_render) {
				array_index = 0;
				long render_sample = 5000000; // 5000000
				double render_light_mult = 1.5; // 1
				int iteration_mult = 1; // 12
				get_fractal(render, fractal_index, is_buddha + is_mandelbulb * 2 + is_mandel4d * 3, z_x, z_y, c_x, c_y, is_mandel ? factor_mandelbrot : factor_julia, power, range, inverse, mult_power, is_mandel, render_width, render_height, iteration >= 50 ? iteration * iteration_mult : iteration, mult_sample * render_sample, min_dist, rot_x, rot_z, raymarch_iterations, phi_shift, fomula_3d, c_z, c_w, z_z, z_w);
				mult_light = light * render_light_mult * (is_buddha ? 0.125 : 1.5);
				for (int y = 0; y < render_height; y++) {
					for (int x = 0; x < render_width; x++) {
						val = render[array_index++];
						if (is_mandelbulb || is_mandel4d) {
							double col = val <= 0 ? 0 : (int)val;
							double ot_div = val - col;
							r = g = b = 0;
							if (ot_div != 0.0) {
								double bright = (1.0 - 1.5 * pow(ot_div, 1.1)) * 255.0;
								vec3 hueo = hue(0.85, 0.45, 0.15, 1.0 + clamp(col, 3.0, 20.0) * 80.0);
								r = clamp(hueo.x * bright, 0, 255);
								g = clamp(hueo.y * bright, 0, 255);
								b = clamp(hueo.z * bright, 0, 255);
							}
						}
						else {
							r = bounce(val * 0.5 * mult_light, 0, 255, 5);
							g = bounce(val * 1.5 * mult_light, 0, 255, 5);
							b = bounce(val * 2.0 * mult_light, 0, 255, 5);
						}
						image.set_pixel(x, y, r, g, b);
					}
				}
				image.save_as("render\\pic" + std::to_string(frame - 1) + ".bmp");
			}
		}

		key = '\b';
		if (__kbhit()) {
			unsigned int start_time = clock();
			key = getch_toupper();
			if (key == '\t') {
				key = getch_toupper();
				repeats_button.push_back(key);
			}
			if (clock() - start_time > 450) {
				repeats_button.clear();
			}
			repeats_button.push_back(key);
			if (key == ' ') {
				repeats_button.clear();
			}
			is_delete_button = 1;
		}
		if (repeats_button.size() > 0) {
			push_button.push_back(repeats_button);
		}
		for (int i = 0; i < repeats_button.size(); i++) {
			char button = repeats_button[i];
			if (type_fractal == 0) {
				if (is_mandelbulb || is_mandel4d) {
					if (button == '2') { move_scale *= 0.5; }
					if (button == '3') { move_scale *= 2.0; }
					if (button == '{') { min_dist *= 0.5; }
					if (button == '}') { min_dist *= 2.0; }
					if (button == 'E') { factor_mandelbrot -= move_scale; }
					if (button == 'Q') { factor_mandelbrot += move_scale; }
					if (button == '6') { rot_z -= factor_mandelbrot * move_scale; }
					if (button == '4') { rot_z += factor_mandelbrot * move_scale; }
					if (button == '8') { rot_x += factor_mandelbrot * move_scale; }
					if (button == '5') { rot_x -= factor_mandelbrot * move_scale; }
					if (button == 'D') { c_x += factor_mandelbrot * move_scale; }
					if (button == 'A') { c_x -= factor_mandelbrot * move_scale; }
					if (button == 'W') { c_y -= factor_mandelbrot * move_scale; }
					if (button == 'S') { c_y += factor_mandelbrot * move_scale; }
					if (button == '(') { phi_speed_mult *= 0.5; }
					if (button == ')') { phi_speed_mult *= 2.0; }
					if (button == '&') { phi_shift += factor_mandelbrot * move_scale * phi_speed_mult; }
					if (button == 'G') { phi_shift = 0.0; }
					if (button == 'I') { fomula_3d = fomula_3d == 0 ? 1 : 0; }
					if (button == '+') { z_z += factor_mandelbrot * move_scale; }
					if (button == '*') { z_z -= factor_mandelbrot * move_scale; }
					if (button == '@') { c_z += factor_mandelbrot * move_scale; }
					if (button == '#') { c_z -= factor_mandelbrot * move_scale; }
					if (button == '?') { z_w += factor_mandelbrot * move_scale; }
					if (button == '>') { z_w -= factor_mandelbrot * move_scale; }
					if (button == '\"') { c_w += factor_mandelbrot * move_scale; }
					if (button == '|') { c_w -= factor_mandelbrot * move_scale; }
				}
				else {
					if (button == 'E') { factor_mandelbrot *= 0.95; }
					if (button == 'Q') { factor_mandelbrot *= 1.05; }
					if (button == 'D') { c_x += factor_mandelbrot * 0.05; }
					if (button == 'A') { c_x -= factor_mandelbrot * 0.05; }
					if (button == 'W') { c_y -= factor_mandelbrot * 0.05; }
					if (button == 'S') { c_y += factor_mandelbrot * 0.05; }
				}
			}
			if (type_fractal == 1) {
				if (button == '6') { c_x += factor_mandelbrot * 0.05; }
				if (button == '4') { c_x -= factor_mandelbrot * 0.05; }
				if (button == '8') { c_y -= factor_mandelbrot * 0.05; }
				if (button == '5') { c_y += factor_mandelbrot * 0.05; }
				if (button == '9') { factor_mandelbrot *= 0.95; }
				if (button == '7') { factor_mandelbrot *= 1.05; }
				if (button == 'D') { z_x += factor_julia * 0.05; }
				if (button == 'A') { z_x -= factor_julia * 0.05; }
				if (button == 'S') { z_y += factor_julia * 0.05; }
				if (button == 'W') { z_y -= factor_julia * 0.05; }
				if (button == 'E') { factor_julia *= 0.95; }
				if (button == 'Q') { factor_julia *= 1.05; }
			}
			if (type_fractal == 2) {
				if (button == 'D') { c_x += factor_mandelbrot * 0.05; }
				if (button == 'A') { c_x -= factor_mandelbrot * 0.05; }
				if (button == 'W') { c_y -= factor_mandelbrot * 0.05; }
				if (button == 'S') { c_y += factor_mandelbrot * 0.05; }
				if (button == 'E') { factor_mandelbrot *= 0.95; }
				if (button == 'Q') { factor_mandelbrot *= 1.05; }
				if (button == '6') { z_x += factor_julia * 0.05; }
				if (button == '4') { z_x -= factor_julia * 0.05; }
				if (button == '8') { z_y -= factor_julia * 0.05; }
				if (button == '5') { z_y += factor_julia * 0.05; }
				if (button == '9') { factor_julia *= 0.95; }
				if (button == '7') { factor_julia *= 1.05; }
			}
			if (button == 'R') { power = (int)(power + 1.0); }
			if (button == 'T') { power = floor(power - 0.001); }
			if (button == 'J') { light *= 0.75; }
			if (button == 'K') { light *= 1.33; }
			if (button == 'M') { scale_mult_power *= 0.5; }
			if (button == ',') { scale_mult_power *= 2.0; }
			if (button == '.') { mult_power -= 0.01 * scale_mult_power; }
			if (button == '/') { mult_power += 0.01 * scale_mult_power; }
			if (button == 'O') { scale_power *= 0.5; }
			if (button == 'P') { scale_power *= 2.0; }
			if (button == '[') { power -= 0.025 * scale_power; }
			if (button == ']') { power += 0.025 * scale_power; }
			if (button == '\'') { inverse -= 0.01 * scale_inverse; }
			if (button == '\\') { inverse += 0.01 * scale_inverse; }
			if (button == 'L') { scale_inverse *= 0.5; }
			if (button == ';') { scale_inverse *= 2.0; }
			if (button == '-') { iteration *= 0.95; }
			if (button == '=') { iteration += 10; }
			if (button == 'Z') { range *= 0.9; }
			if (button == 'X') { range *= 1.1; }
			if (button == 'V') { fractal_index = fractal_index >= (fractal_count - 1) ? 0 : (fractal_index + 1); }
			if (button == 'C') { fractal_index = fractal_index <= 0 ? (fractal_count - 1) : (fractal_index - 1); }
			if (button == '`') { power = 2.0, range = 8.0, inverse = 0.0, mult_power = 1.0, scale_mult_power = 1.0, scale_inverse = 1.0, scale_power = 1.0, factor_mandelbrot = 1.5, factor_julia = 1.5, type_fractal = 0; is_mandel = 1; phi_shift = 0.0; }
			if (button == '3') { mult_sample *= 0.5; }
			if (button == '2') { mult_sample *= 2.0; }
			if (button == 'F') { inverse = inverse == 1.0 ? 0.0 : 1.0; }
			if (button == 'N') { is_mandel4d = 0; is_mandelbulb = is_mandelbulb == 0 ? 1 : 0; is_buddha = 0; }
			if (button == 'Y') { is_mandelbulb = 0; is_mandel4d = is_mandel4d == 0 ? 1 : 0; is_buddha = 0; }
			if (button == 'B') { is_buddha = is_buddha == 0 ? 1 : 0; is_mandelbulb = 0; }
			if (button == '0') { is_render = is_render == 0 ? 1 : 0; }
			if (button == ' ') { repeats_button.clear(); break; }
			if (button == '!') { is_ansi = is_ansi == 0 ? 1 : 0; clear_screen(); }
			if (button == 'U') { type_fractal = (type_fractal + 2) % 3; is_mandel = (type_fractal == 0); }
			if (button == 'H') {
				if (type_fractal == 0 || type_fractal == 2) { c_x = c_y = c_z = c_w = rot_x = rot_z = 0.0; factor_mandelbrot = 1.5; }
				if (type_fractal == 1 || type_fractal == 2) { z_x = z_y = z_z = z_w = rot_x = rot_z = 0.0; factor_julia = 1.5; }
			}
			if (button == '_' && frame > 30) {
				std::ofstream fout("button.txt");
				for (const auto& x: push_button) {
					for (const char ch : x) {
						if (ch != '_') {
							fout << ch;
						}
					}
					fout << '\n';
				}
				fout.close();
			}
			if (button == '~' && frame <= 30) {
				frame = 1;
				std::ifstream fin("button.txt");
				std::string line;
				while (std::getline(fin, line)) {
					sequence_button.emplace_back(line.begin(), line.end());
				}
				fin.close();
			}
		}
		if (repeats_button.size() > 0 && is_delete_button) {
			repeats_button.pop_back();
			is_delete_button = 0;
		}
		if (sequence_button.size() > 0) {
			repeats_button = sequence_button[0];
			sequence_button.erase(sequence_button.begin());
			is_delete_button = 1;
		}
	}
}