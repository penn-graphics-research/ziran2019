#ifndef PRETTY_PRINTING
#define PRETTY_PRINTING

#include <iostream>
#include <chrono>
#include <string>

namespace ZIRAN {

/**
 * Enums: Color codes
 *
 * BLACK - Black
 * BLUE - Blue
 * GREEN - Green
 * CYAN - Cyan
 * RED - Red
 * MAGENTA - Magenta / purple
 * BROWN - Brown / dark yellow
 * GREY - Grey / dark white
 * DARKGREY - Dark grey / light black
 * LIGHTBLUE - Light blue
 * LIGHTGREEN - Light green
 * LIGHTCYAN - Light cyan
 * LIGHTRED - Light red
 * LIGHTMAGENTA - Light magenta / light purple
 * YELLOW - Yellow (bright)
 * WHITE - White (bright)
 * RESET - Reset to default
 */
enum class Color {
    BLACK,
    BLUE,
    GREEN,
    CYAN,
    RED,
    MAGENTA,
    BROWN,
    GREY,
    DARKGREY,
    LIGHTBLUE,
    LIGHTGREEN,
    LIGHTCYAN,
    LIGHTRED,
    LIGHTMAGENTA,
    YELLOW,
    WHITE,
    RESET,
    NOCHANGE
};

/**
 * Consts: ANSI escape strings
 *
 * ANSI_CLS                - Clears screen
 * ANSI_CONSOLE_TITLE_PRE  - Prefix for changing the window title, print the window title in between
 * ANSI_CONSOLE_TITLE_POST - Suffix for changing the window title, print the window title in between
 * ANSI_ATTRIBUTE_RESET    - Resets all attributes
 * ANSI_CURSOR_HIDE        - Hides the cursor
 * ANSI_CURSOR_SHOW        - Shows the cursor
 * ANSI_CURSOR_HOME        - Moves the cursor home (0,0)
 * ANSI_BLACK              - Black
 * ANSI_RED                - Red
 * ANSI_GREEN              - Green
 * ANSI_BROWN              - Brown / dark yellow
 * ANSI_BLUE               - Blue
 * ANSI_MAGENTA            - Magenta / purple
 * ANSI_CYAN               - Cyan
 * ANSI_GREY               - Grey / dark white
 * ANSI_DARKGREY           - Dark grey / light black
 * ANSI_LIGHTRED           - Light red
 * ANSI_LIGHTGREEN         - Light green
 * ANSI_YELLOW             - Yellow (bright)
 * ANSI_LIGHTBLUE          - Light blue
 * ANSI_LIGHTMAGENTA       - Light magenta / light purple
 * ANSI_LIGHTCYAN          - Light cyan
 * ANSI_WHITE              - White (bright)
 * ANSI_BACKGROUND_BLACK   - Black background
 * ANSI_BACKGROUND_RED     - Red background
 * ANSI_BACKGROUND_GREEN   - Green background
 * ANSI_BACKGROUND_YELLOW  - Yellow background
 * ANSI_BACKGROUND_BLUE    - Blue background
 * ANSI_BACKGROUND_MAGENTA - Magenta / purple background
 * ANSI_BACKGROUND_CYAN    - Cyan background
 * ANSI_BACKGROUND_WHITE   - White background
 */
// Remaining colors not supported as background colors
std::ostream& operator<<(std::ostream& os, const Color& c);

std::string stripColors(const std::string& s);

struct SetColor {
    Color fg, bg;

    SetColor(Color fg, Color bg = Color::NOCHANGE);

    std::ostream& operator()(std::ostream& os) const;

    std::ostream& setForegroundColor(std::ostream& os) const;

    std::ostream& setBackgroundColor(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const SetColor& p);

// Needed for the windows console to support ansi color codes
int initializeColors();

std::string setConsoleTitle(const std::string& title);

std::string setCursorVisiblity(bool visible);

void showCursor(bool visible);

void hideCursor(bool visible);

std::string moveCursor(int line, int column);

std::string moveCursorRelative(int lines_down, int columns_right);

std::string clearScreen();

void getTerminalSize(int& lines, int& columns);

struct ProgressBar {
    double progress;
    int barWidth;
    std::string tail;
    ProgressBar(double progress, const std::string& tail = "", int barWidth = 73);
    std::ostream& operator()(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const ProgressBar& p);

struct PrettyPrintDuration {
    std::chrono::duration<double> elapsed_seconds;
    PrettyPrintDuration(std::chrono::duration<double> elapsed_seconds);

    std::ostream& operator()(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const PrettyPrintDuration& p);
} // namespace ZIRAN
#endif /* ifndef PRETTY_PRINTING */
