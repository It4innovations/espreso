/**
 * @file
 *  This file is part of UTILS
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2014-2015, Technische Universitaet Muenchen.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from this
 *     software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef UTILS_PROGRESS_H
#define UTILS_PROGRESS_H

#include "utils/env.h"
#include "utils/logger.h"
#include "utils/stringutils.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>
#include <sys/ioctl.h>

namespace utils
{

// TODO replace this with constexpr in C++11
#define ROTATION_IND "-\\|/"

class Progress
{
private:
	enum OutputType
	{
		DISABLED,
		TTY
		// TODO Support file output
	};

private:
	/** The output stream we use */
	std::ostream *m_output;

	OutputType m_type;

	/** Total number of updates */
	unsigned long m_total;

	/** Current update status */
	unsigned long m_current;

	/** Size of the progress bar */
	unsigned long m_barSize;

	/** Rotation indicator position */
	unsigned char m_rotPosition;

	/** TTY handle (if used) */
	std::ofstream m_tty;

public:
	Progress(unsigned long total = 100)
		: m_current(0), m_barSize(80) /* Default size */, m_rotPosition(0)
	{
		std::string envOutput = Env::get<std::string>("UTILS_PROGRESS_OUTPUT", "STDERR");

		StringUtils::toUpper(envOutput);

		if (envOutput == "STDOUT") {
			m_output = &std::cout;
			m_type = TTY;

			setSize(isatty(fileno(stdout)));
		} else if (envOutput == "STDERR") {
			m_output = &std::cerr;
			m_type = TTY;

			setSize(isatty(fileno(stderr)));
		} else if (envOutput == "TTY") {
			m_tty.open("/dev/tty"); // try unix
			if (!m_tty)
				m_tty.open("CON:"); // try windows

			if (m_tty) {
				m_output = &m_tty;
				m_type = TTY;
				setSize();
			} else {
				logWarning() << "Could not open terminal. Disabling progress bar.";
				m_type = DISABLED;
			}
		} else {
			m_type = DISABLED;
		}
	}

	/**
	 * Set a new total value
	 * Does not update the progress bar
	 */
	void setTotal(unsigned long total)
	{
		m_total = total;
	}

	/**
	 * Set the current value of the progress bar without updating the screen
	 */
	void set(unsigned long current)
	{
		m_current = std::min(current, m_total);
	}

	/**
	 * Update the progress bar
	 *
	 * If <code>current >= total</code> a newline is printed and it is no longer possible
	 * to update screen
	 */
	void update(unsigned long current)
	{
		set(current);

		if (m_type == DISABLED)
			return;

		// Calculuate the ratio of complete-to-incomplete.
		const float ratio = m_current/static_cast<float>(m_total);

		// Show the percentage complete
		(*m_output) << std::setw(3) << (int)(ratio*100) << "% [";

		// real width (without additional chars)
		const unsigned long realSize = m_barSize - 9;

		const unsigned long comChars = realSize * ratio;

		// Show the load bar
		for (unsigned int i = 0; i < comChars; i++)
			(*m_output) << '=';

		for (unsigned int i = comChars; i < realSize; i++)
			(*m_output) << ' ';

		(*m_output) << "] ";

		// Print rotation indicator
		(*m_output) << ROTATION_IND[m_rotPosition];
		m_rotPosition = (m_rotPosition+1) % 4;

		// go to the beginning of the line
		(*m_output) << '\r' << std::flush;
	}

	/**
	 * Updates the progress bar but does not change the current value
	 */
	void update()
	{
		update(m_current);
	}

	/**
	 * Updates the progress bar and increments it by one
	 */
	void increment()
	{
		update(m_current + 1);
	}

	/**
	 * Removes the progress bar from the output
	 */
	void clear()
	{
		if (m_type == DISABLED)
			return;

		for (unsigned int i = 0; i < m_barSize; i++)
			(*m_output) << ' ';
		(*m_output) << '\r' << std::flush;
	}

private:
	/**
	 * Sets progress bar size according to the terminal size
	 */ 
	void setSize(bool automatic = true)
	{
		// Check if size is set in env
		unsigned long size = Env::get<unsigned long>("UTILS_PROGRESS_SIZE", 0);

		if (size > 0) {
			m_barSize = size;
			return;
		}

		if (!automatic)
			// No automatic size detection (e.g. for stdout)
			return;

		// Try to get terminal size
#ifdef TIOCGSIZE
		struct ttysize ts;
		ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
		m_barSize = ts.ts_cols;
#elif defined(TIOCGWINSZ)
		struct winsize ts;
		ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
		m_barSize = ts.ws_col;
#else
		logWarning() << "Could not get terminal size, using default";
#endif
	}
};

}

#endif // UTILS_PROGRESS_H
