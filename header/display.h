/*
 * display.h
 *
 *  Created on: Nov 4, 2016
 *      Author: bryan
 */

#ifndef DISPLAY_H_
#define DISPLAY_H_

#include <string>
#include <SDL2/SDL.h>

class Display
{
public:
	Display(int width, int height, const std::string& title);

	void Clear(float r, float g, float b, float a);
	void update();
	bool IsClosed();

	virtual ~Display();
protected:
private:
	void operator=(const Display& display) {}
	Display(const Display& display) {}

	SDL_Window* m_window;
	SDL_GLContext m_glContext;
	bool m_isClosed;
};


#endif /* DISPLAY_H_ */
