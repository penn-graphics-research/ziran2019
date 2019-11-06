#include <Ziran/CS/Util/PrettyPrinting.h>

#include <iomanip>
#include <regex>
#ifndef _WIN32
#include <sys/ioctl.h>
#include <unistd.h>
#else
#include <WinSock2.h>
#include <Windows.h>
#include <io.h>
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define isatty _isatty
#endif

namespace ZIRAN {

// Remaining colors not supported as background colors
std::ostream& operator<<(std::ostream& os, const Color& c)
{
    const char* ANSI_ATTRIBUTE_RESET = "\033[0m";
    const char* ANSI_BLACK = "\033[22;30m";
    const char* ANSI_RED = "\033[22;31m";
    const char* ANSI_GREEN = "\033[22;32m";
    const char* ANSI_BROWN = "\033[22;33m";
    const char* ANSI_BLUE = "\033[22;34m";
    const char* ANSI_MAGENTA = "\033[22;35m";
    const char* ANSI_CYAN = "\033[22;36m";
    const char* ANSI_GREY = "\033[22;37m";
    const char* ANSI_DARKGREY = "\033[01;30m";
    const char* ANSI_LIGHTRED = "\033[01;31m";
    const char* ANSI_LIGHTGREEN = "\033[01;32m";
    const char* ANSI_YELLOW = "\033[01;33m";
    const char* ANSI_LIGHTBLUE = "\033[01;34m";
    const char* ANSI_LIGHTMAGENTA = "\033[01;35m";
    const char* ANSI_LIGHTCYAN = "\033[01;36m";
    const char* ANSI_WHITE = "\033[01;37m";
    switch (c) {
    case Color::BLACK:
        os << ANSI_BLACK;
        break;
    case Color::BLUE:
        os << ANSI_BLUE;
        break;
    case Color::GREEN:
        os << ANSI_GREEN;
        break;
    case Color::CYAN:
        os << ANSI_CYAN;
        break;
    case Color::RED:
        os << ANSI_RED;
        break;
    case Color::MAGENTA:
        os << ANSI_MAGENTA;
        break;
    case Color::BROWN:
        os << ANSI_BROWN;
        break;
    case Color::GREY:
        os << ANSI_GREY;
        break;
    case Color::DARKGREY:
        os << ANSI_DARKGREY;
        break;
    case Color::LIGHTBLUE:
        os << ANSI_LIGHTBLUE;
        break;
    case Color::LIGHTGREEN:
        os << ANSI_LIGHTGREEN;
        break;
    case Color::LIGHTCYAN:
        os << ANSI_LIGHTCYAN;
        break;
    case Color::LIGHTRED:
        os << ANSI_LIGHTRED;
        break;
    case Color::LIGHTMAGENTA:
        os << ANSI_LIGHTMAGENTA;
        break;
    case Color::YELLOW:
        os << ANSI_YELLOW;
        break;
    case Color::WHITE:
        os << ANSI_WHITE;
        break;
    case Color::RESET:
        os << ANSI_ATTRIBUTE_RESET;
        break;
    case Color::NOCHANGE:
        break;
    }
    return os;
}

std::string stripColors(const std::string& s)
{
    static const std::regex r("[\033\007][[()#;?]*(?:[0-9]{1,4}(?:;[0-9]{0,4})*)?[0-9A-ORZcf-nqry=><]");
    return std::regex_replace(s, r, "");
}

SetColor::SetColor(Color fg, Color bg)
    : fg(fg)
    , bg(bg)
{
}

std::ostream& SetColor::operator()(std::ostream& os) const
{
    const char* ANSI_ATTRIBUTE_RESET = "\033[0m";
    if (bg == Color::RESET) {
        os << ANSI_ATTRIBUTE_RESET;
        setForegroundColor(os);
    }
    else {
        setForegroundColor(os);
        setBackgroundColor(os);
    }
    return os;
}

std::ostream& SetColor::setForegroundColor(std::ostream& os) const
{
    return os << fg;
}

std::ostream& SetColor::setBackgroundColor(std::ostream& os) const
{
    const char* ANSI_BACKGROUND_BLACK = "\033[40m";
    const char* ANSI_BACKGROUND_RED = "\033[41m";
    const char* ANSI_BACKGROUND_GREEN = "\033[42m";
    const char* ANSI_BACKGROUND_YELLOW = "\033[43m";
    const char* ANSI_BACKGROUND_BLUE = "\033[44m";
    const char* ANSI_BACKGROUND_MAGENTA = "\033[45m";
    const char* ANSI_BACKGROUND_CYAN = "\033[46m";
    const char* ANSI_BACKGROUND_WHITE = "\033[47m";
    const char* ANSI_ATTRIBUTE_RESET = "\033[0m";
    switch (bg) {
    case Color::BLACK:
        os << ANSI_BACKGROUND_BLACK;
        break;
    case Color::BLUE:
        os << ANSI_BACKGROUND_BLUE;
        break;
    case Color::GREEN:
        os << ANSI_BACKGROUND_GREEN;
        break;
    case Color::CYAN:
        os << ANSI_BACKGROUND_CYAN;
        break;
    case Color::RED:
        os << ANSI_BACKGROUND_RED;
        break;
    case Color::MAGENTA:
        os << ANSI_BACKGROUND_MAGENTA;
        break;
    case Color::BROWN:
        os << ANSI_BACKGROUND_YELLOW;
        break;
    case Color::GREY:
        os << ANSI_BACKGROUND_WHITE;
        break;
    case Color::RESET:
        os << ANSI_ATTRIBUTE_RESET;
        break;
    default:
        break;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const SetColor& p)
{
    return p(os);
}

// Needed for the windows console to support ansi color codes
int initializeColors()
{
#ifdef _WIN32
    // Set output mode to handle virtual terminal sequences
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hOut == INVALID_HANDLE_VALUE) {
        return GetLastError();
    }

    DWORD dwMode = 0;
    if (!GetConsoleMode(hOut, &dwMode)) {
        return GetLastError();
    }

    dwMode |= 0x0004; //ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    if (!SetConsoleMode(hOut, dwMode)) {
        return GetLastError();
    }
#endif
    return 0;
}

std::string setConsoleTitle(const std::string& title)
{
    const char ANSI_CONSOLE_TITLE_PRE[] = "\033]0;";
    const char ANSI_CONSOLE_TITLE_POST[] = "\007";
    return ANSI_CONSOLE_TITLE_PRE + title + ANSI_CONSOLE_TITLE_POST;
}

std::string setCursorVisiblity(bool visible)
{
    const char ANSI_CURSOR_HIDE[] = "\033[?25l";
    const char ANSI_CURSOR_SHOW[] = "\033[?25h";
    return (visible ? ANSI_CURSOR_SHOW : ANSI_CURSOR_HIDE);
}

void showCursor(bool visible)
{
    setCursorVisiblity(true);
}

void hideCursor(bool visible)
{
    setCursorVisiblity(false);
}

std::string moveCursor(int line, int column)
{
    return "\033[" + std::to_string(line) + ';' + std::to_string(column) + 'H';
}

std::string moveCursorRelative(int lines_down, int columns_right)
{
    return "\033[" + std::to_string(abs(lines_down)) + (lines_down < 0 ? 'A' : 'B')
        + "\033[" + std::to_string(abs(columns_right)) + (columns_right < 0 ? 'D' : 'C');
}

std::string clearScreen()
{
    return std::string("\033[2J\033[3J\033[H");
}

void getTerminalSize(int& lines, int& columns)
{
#ifdef TIOCGSIZE
    struct ttysize ts;
    ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
    lines = ts.ts_lines;
    columns = ts.ts_cols;
#elif defined TIOCGWINSZ
    struct winsize ts;
    ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
    lines = ts.ws_row;
    columns = ts.ws_col;
#endif
}

ProgressBar::ProgressBar(double progress, const std::string& tail, int barWidth)
    : progress(progress)
    , barWidth(barWidth)
    , tail(tail)
{
}

std::ostream& ProgressBar::operator()(std::ostream& os) const
{
    os << "[";
    os << Color::GREEN;
    int pos = int(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            os << "=";
        else if (i == pos)
            os << ">";
        else
            os << " ";
    }
    os << Color::RESET;
    os << "] ";
    os << std::setprecision(3) << std::setw(6) << std::fixed << progress * 100.0 << " " << tail;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ProgressBar& p)
{
    return p(os);
}

PrettyPrintDuration::PrettyPrintDuration(std::chrono::duration<double> elapsed_seconds)
    : elapsed_seconds(elapsed_seconds)
{
}
std::ostream& PrettyPrintDuration::operator()(std::ostream& os) const
{
    using namespace std::chrono;
    auto dur = elapsed_seconds;
    int num_days = int(dur / hours(24));
    if (num_days)
        os << Color::RED << num_days << "d ";
    dur -= num_days * hours(24);
    int num_hours = int(dur / hours(1));
    if (num_hours && !num_days)
        os << Color::BROWN;
    if (num_hours)
        os << num_hours << "h ";
    dur -= num_hours * hours(1);

    int num_minutes = int(dur / minutes(1));
    if (!num_hours)
        os << Color::GREEN;
    if (num_minutes)
        os << num_minutes << "m ";
    dur -= num_minutes * minutes(1);
    os << dur.count() << "s" << Color::RESET;
    return os;
}

std::ostream& operator<<(std::ostream& os, const PrettyPrintDuration& p)
{
    return p(os);
}
} // namespace ZIRAN
