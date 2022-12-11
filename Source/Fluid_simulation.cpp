#include <iostream>
#include <SFML/Graphics.hpp>
#include <windows.h>

using namespace sf;



constexpr auto ITER = 6;		// кількість ітерацій для метода solve()
constexpr auto N = 256;			// кількість клітинок у вікні 256х256
constexpr auto SCALE = 2;		// збільшення клітинок

int INDEX(int x, int y) {		// метод для повернення 1-вимірної координати в 2-вимірній сітці
	x = x > N - 1 ?	 N - 1 : x;
	y = y > N - 1 ? N - 1 : y;
	
	x = x < 0 ? 0 : x;
	y = y < 0 ? 0 : y;
	
	return (x + y * N);
}

void set_bound(int b, float* x)		//метод для встановлення обмежувачів за які рідина не витікатиме 
{

	for (int i = 1; i < N - 1; i++) {
		x[INDEX(i, 0)] = b == 2 ? -x[INDEX(i, 1)] : x[INDEX(i, 1)];
		x[INDEX(i, N - 1)] = b == 2 ? -x[INDEX(i, N - 2)] : x[INDEX(i, N - 2)];
	}


	for (int i = 1; i < N - 1; i++) {
		x[INDEX(0, i)] = b == 1 ? -x[INDEX(1, i)] : x[INDEX(1, i)];
		x[INDEX(N - 1, i)] = b == 1 ? -x[INDEX(N - 2, i)] : x[INDEX(N - 2, i)];
	}


	x[INDEX(0, 0)] = 0.5f * (x[INDEX(1, 0)] + x[INDEX(0, 1)]);
	x[INDEX(0, N - 1)] = 0.5f * (x[INDEX(1, N - 1)] + x[INDEX(0, N - 2)]);
	x[INDEX(N - 1, 0)] = 0.5f * (x[INDEX(N - 1, 1)] + x[INDEX(N - 2, 0)]);
	x[INDEX(N - 1, N - 1)] = 0.5f * (x[INDEX(N - 2, N - 1)] + x[INDEX(N - 1, N - 2)]);
}

void solve(int b, float* x, float* x0, float a, float c) {			// метод для визначення приблизних рішень лінійного рівняння чим більше ITER, тим більша точність
	float cRecip = 1.f / c;
	for (int k = 0; k < ITER; k++) {
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				x[INDEX(i, j)] = 
					(x0[INDEX(i, j)] + a * (x[INDEX(i + 1, j)] + x[INDEX(i - 1, j)] 
						+ x[INDEX(i, j + 1)] + x[INDEX(i, j - 1)])) * cRecip;
			}
		}

		set_bound(b, x);
	}

}

void diffuse(int b, float* x, float* x0, float diff, float dt)		// метод для дифузії по будь-якій величині яка має початкове значення х0, та наступне х
{
	float a = dt * diff * (N - 2) * (N - 2);
	solve(b, x, x0, a, 1 + 4 * a);
}

void project(float* velocX, float* velocY, float* p, float* div)	// метод для того щоб привести швидкості в умову нестисливості,
{																	// швидкості не з'явилися з нічого і не зникли в нікуди (div(V) = 0)

	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			div[INDEX(i, j)] = (-0.5f * (
				velocX[INDEX(i + 1, j)]
				- velocX[INDEX(i - 1, j)]
				+ velocY[INDEX(i, j + 1)]
				- velocY[INDEX(i, j - 1)])) / N;
			p[INDEX(i, j)] = 0;
		}
	}

	set_bound(0, div);
	set_bound(0, p);
	solve(0, p, div, 1, 6);


	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			velocX[INDEX(i, j)] -= 0.5f * (p[INDEX(i + 1, j)] - p[INDEX(i - 1, j)]) * N;
			velocY[INDEX(i, j)] -= 0.5f * (p[INDEX(i, j + 1)] - p[INDEX(i, j - 1)]) * N;

		}
	}

	set_bound(1, velocX);
	set_bound(2, velocY);
}

void advect(int b, float* d, float* d0, float* velocX, float* velocY, float dt)		// метод для переміщення швидкостей та "фарби" згідно дифузії
{
	float i0, i1, j0, j1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);


	float s0, s1, t0, t1;
	float tmp1, tmp2, x, y;

	float Nfloat = N - 2;
	float ifloat, jfloat;
	int i, j;

	for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
		for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
			tmp1 = dtx * velocX[INDEX(i, j)];
			tmp2 = dty * velocY[INDEX(i, j)];

			x = ifloat - tmp1;
			y = jfloat - tmp2;


			if (x < 0.5f) x = 0.5f;
			if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
			i0 = floorf(x);
			i1 = i0 + 1.0f;
			if (y < 0.5f) y = 0.5f;
			if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
			j0 = floorf(y);
			j1 = j0 + 1.0f;

			s1 = x - i0;
			s0 = 1.0f - s1;
			t1 = y - j0;
			t0 = 1.0f - t1;

			int i0_int = (int)i0;
			int i1_int = (int)i1;
			int j0_int = (int)j0;
			int j1_int = (int)j1;

			d[INDEX(i, j)] =
				s0 * (t0 * d0[INDEX(i0_int, j0_int)] + t1 * d0[INDEX(i0_int, j1_int)])
			   +s1 * (t0 * d0[INDEX(i1_int, j0_int)] + t1 * d0[INDEX(i1_int, j1_int)]);
		}
	}

	set_bound(b, d);
}







class FluidSquare {													// клас рідини
public: FluidSquare(float diff, float visc, float dt) {				// конструктор класу

		this->diffusion = diff;
		this->viscosity = visc;
		this->dt = dt;


		this->prev_dens = new float[N * N];
		this->density = new float[N * N];

		this->Vx = new float[N * N];
		this->Vy = new float[N * N];

		this->Vx0 = new float[N * N];
		this->Vy0 = new float[N * N];
		for (int i = 0; i < N * N; i++) {
			prev_dens[i] = 0;
			density[i] = 0;
			Vx[i] = 0;
			Vy[i] = 0;
			Vx0[i] = 0;
			Vy0[i] = 0;
		}
	}

	~FluidSquare() {												// деструктор класу
		delete[] prev_dens;
		delete[] density;
		delete[] Vx;
		delete[] Vy;
		delete[] Vx0;
		delete[] Vy0;
	}

	float dt;			// крок часу
	float diffusion;	// дифузія
	float viscosity;	// в'язкість

	float* prev_dens;	// попередня в'язкість
	float* density;		// теперішня в'язкість

	float* Vx;			// теперішня швидкість по іксу
	float* Vy;			// теперішня швидкість по ігреку

	float* Vx0;			// попередня швидкість по іксу
	float* Vy0;			// попередня швидкість по ігреку

	void FluidStep()	// метод "кроку" рідини
	{

		diffuse(1, Vx0, Vx, viscosity, dt);		// робимо дифузію х складових швидкостей
		diffuse(2, Vy0, Vy, viscosity, dt);		// робимо дифузію у складових швидкостей

		project(Vx0, Vy0, Vx, Vy);				// приводимо все до нестисливого вигляду

		advect(1, Vx, Vx0, Vx0, Vy0, dt);		// переміщуємо х складові швидкості згідно дифузії
		advect(2, Vy, Vy0, Vx0, Vy0, dt);		// переміщуємо у складові швидкості згідно дифузії

		project(Vx, Vy, Vx0, Vy0);				// приводимо нові значення швидкосией до нестисливого вигляду

		diffuse(0, prev_dens, density, diffusion, dt);		// робимо дифузію густини у рідині
		advect(0, density, prev_dens, Vx, Vy, dt);			// переміщуємо густини по рідині згідно дифузії
	}

	void AddDensity(int x, int y, float amount) {						// метод для додавання "фарби" у точку (х, у)
		int index = INDEX(x, y);
		this->density[index] += amount;
	}

	void AddVelocity(int x, int y, float amountX, float amountY) {		// метод для додавання х та у складової швидкості amountX та amountY у точку (х, у) відповідно
		int index = INDEX(x, y);
		this->Vx[index] += amountX;
		this->Vy[index] += amountY;
	}


};









int main() {

	RenderWindow window(VideoMode(N*SCALE, N*SCALE), "Fluid Simulation");	// створюємо вікно 
	
	RectangleShape* cells_array = new RectangleShape[N*N];					// задаємо масив клітинок
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			
			cells_array[INDEX(i, j)].setSize(Vector2f(SCALE, SCALE));
			cells_array[INDEX(i, j)].setPosition(Vector2f(i*SCALE, j*SCALE));
			
		}
	}
	
	
	

	


	bool is_vel = false, is_V_pressed = false;

	FluidSquare fluid(0.000001f, 0.00007f, 0.01f);				// створюємо об'єкт класу рідини з заданими параметрами дифузії швидкості, в'язеості та кроку часу
	

	
	while (window.isOpen())										// головний цикл програми. Виконується, поки вікно відкрито
	{
		
		Event event;											// обробка подій у вікні
		while (window.pollEvent(event))
		{
			
			if (event.type == Event::Closed)					// якщо натиснуто "хрестик" - закрити вікно
				
				window.close();
		}
		
		
		
		

		if (Keyboard::isKeyPressed(Keyboard::C)){				// якщо натиснуто "С" - очистити рідину від "фарби"
			for (int i = 0; i < N*N; i++)
			{
				fluid.prev_dens[i] = 0;
				fluid.density[i] = 0;
				fluid.Vx0[i] = 0;
				fluid.Vx[i] = 0;
				fluid.Vy0[i] = 0;
				fluid.Vy[i] = 0;

			}
		}
		
		if (is_V_pressed) {
			if (Event::KeyReleased) {
				if (event.key.code == Keyboard::V) {
					is_vel = is_vel ? false : true;
					is_V_pressed = false;
				}
			}
		}
		if (Keyboard::isKeyPressed(Keyboard::V)) {
			is_V_pressed = true;
		}


		int amount = 400;
		if (Mouse::isButtonPressed(Mouse::Left)) {				// якщо натиснуто ліву кнопку миші - додати "фарбу" у позицію курсора та область поряд з нею
																// та додати швидкість згідно напрямку руху курсора
																
			
			int x = Mouse::getPosition(window).x;
			int y = Mouse::getPosition(window).y;
			
			fluid.AddDensity(x / SCALE, y / SCALE, 4 * amount);
			
			fluid.AddDensity(x / SCALE -1, y / SCALE, 2 * amount);
			fluid.AddDensity(x / SCALE +1, y / SCALE, 2 * amount);
			fluid.AddDensity(x / SCALE, y / SCALE -1, 2 * amount);
			fluid.AddDensity(x / SCALE, y / SCALE +1, 2 * amount);
			
			fluid.AddDensity(x / SCALE - 1, y / SCALE - 1, 2 * amount);
			fluid.AddDensity(x / SCALE + 1, y / SCALE + 1, 2 * amount);
			fluid.AddDensity(x / SCALE + 1, y / SCALE - 1, 2 * amount);
			fluid.AddDensity(x / SCALE - 1, y / SCALE + 1, 2 * amount);
			fluid.AddDensity(x / SCALE - 2, y / SCALE, 2 * amount);
			fluid.AddDensity(x / SCALE + 2, y / SCALE, 2 * amount);
			fluid.AddDensity(x / SCALE, y / SCALE - 2, 2 * amount);
			fluid.AddDensity(x / SCALE, y / SCALE + 2, 2 * amount);
			
			fluid.AddDensity(x / SCALE - 3, y / SCALE, 1 * amount);
			fluid.AddDensity(x / SCALE + 3, y / SCALE, 1 * amount);
			fluid.AddDensity(x / SCALE, y / SCALE - 3, 1 * amount);
			fluid.AddDensity(x / SCALE, y / SCALE + 3, 1 * amount);
			fluid.AddDensity(x / SCALE + 2, y / SCALE + 1, 1 * amount);
			fluid.AddDensity(x / SCALE + 2, y / SCALE - 1, 1 * amount);
			fluid.AddDensity(x / SCALE - 2, y / SCALE + 1, 1 * amount);
			fluid.AddDensity(x / SCALE - 2, y / SCALE - 1, 1 * amount);
			fluid.AddDensity(x / SCALE - 1, y / SCALE + 2, 1 * amount);
			fluid.AddDensity(x / SCALE + 1, y / SCALE + 2, 1 * amount);
			fluid.AddDensity(x / SCALE - 1, y / SCALE - 2, 1 * amount);
			fluid.AddDensity(x / SCALE + 1, y / SCALE - 2, 1 * amount);
			fluid.AddDensity(x / SCALE - 2, y / SCALE - 2, 1 * amount);
			fluid.AddDensity(x / SCALE + 2, y / SCALE + 2, 1 * amount);
			fluid.AddDensity(x / SCALE + 2, y / SCALE - 2, 1 * amount);
			fluid.AddDensity(x / SCALE - 2, y / SCALE + 2, 1 * amount);
			
			Sleep(1);																//затримка в 1 мс для визначення переміщення курсору
			fluid.AddVelocity(x / SCALE, y / SCALE, 500*(Mouse::getPosition(window).x  - x), 500*(Mouse::getPosition(window).y  - y));
				
		}




		if (Mouse::isButtonPressed(Mouse::Right)) {									// якщо натиснуто праву кнопку миші - додати швидкість згідно напрямку руху курсора
			int x = Mouse::getPosition(window).x;
			int y = Mouse::getPosition(window).y;
			Sleep(1);																// затримка в 1 мс для визначення переміщення курсору
			fluid.AddVelocity(x / SCALE, y / SCALE, 500 * (Mouse::getPosition(window).x - x), 500 * (Mouse::getPosition(window).y - y));
		}


		fluid.FluidStep();															// крок рідини
		


		window.clear();																// очищення екрану від минулого кадру


		for (int i = 0; i < N*N; i++) {												// фарбування клітинок на екрані за густиною клітинок рідини
			float brightness = fluid.density[i];
			if (brightness > 255) brightness = 255;
			cells_array[i].setFillColor(Color(int(brightness), 102, 255, int(brightness)));
			window.draw(cells_array[i]);											
		}
		
		if (is_vel) {
			for (int j = 0; j < N / 4; j++) {
				for (int i = 0; i < N / 4; i++) {
					float midx =
						(fluid.Vx[i * 4 + j * 4 * N] + fluid.Vx[i * 4 + 1 + j * 4 * N] + fluid.Vx[i * 4 + 2 + j * 4 * N] + fluid.Vx[i * 4 + 3 + j * 4 * N] +
							fluid.Vx[i * 4 + j * 4 * N + 1] + fluid.Vx[i * 4 + 1 + j * 4 * N + 1] + fluid.Vx[i * 4 + 2 + j * 4 * N + 1] + fluid.Vx[i * 4 + 3 + j * 4 * N + 1] +
							fluid.Vx[i * 4 + j * 4 * N + 2] + fluid.Vx[i * 4 + 1 + j * 4 * N + 2] + fluid.Vx[i * 4 + 2 + j * 4 * N + 2] + fluid.Vx[i * 4 + 3 + j * 4 * N + 2] +
							fluid.Vx[i * 4 + j * 4 * N + 3] + fluid.Vx[i * 4 + 1 + j * 4 * N + 3] + fluid.Vx[i * 4 + 2 + j * 4 * N + 3] + fluid.Vx[i * 4 + 3 + j * 4 * N + 3]) / 4;
					float midy =
						(fluid.Vy[i * 4 + j * 4 * N] + fluid.Vy[i * 4 + 1 + j * 4 * N] + fluid.Vy[i * 4 + 2 + j * 4 * N] + fluid.Vy[i * 4 + 3 + j * 4 * N] +
							fluid.Vy[i * 4 + j * 4 * N + 1] + fluid.Vy[i * 4 + 1 + j * 4 * N + 1] + fluid.Vy[i * 4 + 2 + j * 4 * N + 1] + fluid.Vy[i * 4 + 3 + j * 4 * N + 1] +
							fluid.Vy[i * 4 + j * 4 * N + 2] + fluid.Vy[i * 4 + 1 + j * 4 * N + 2] + fluid.Vy[i * 4 + 2 + j * 4 * N + 2] + fluid.Vy[i * 4 + 3 + j * 4 * N + 2] +
							fluid.Vy[i * 4 + j * 4 * N + 3] + fluid.Vy[i * 4 + 1 + j * 4 * N + 3] + fluid.Vy[i * 4 + 2 + j * 4 * N + 3] + fluid.Vy[i * 4 + 3 + j * 4 * N + 3]) / 4;
					Vertex line[] =
					{
						Vertex(Vector2f(SCALE * 4 * i + 4, SCALE * 4 * j + 4)),
						Vertex(Vector2f(SCALE * 4 * (midx + i) + 4, SCALE * 4 * (midy + j) + 4))
					};
					window.draw(line, 2, sf::Lines);
				}
			}
		}
		
		window.display();															// відображення намальованих клітинок на екрані
	}

	delete[] cells_array;															
	return 0;
}
