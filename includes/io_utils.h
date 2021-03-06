#ifndef _IO_UTILS_H_
#define _IO_UTILS_H_

/* If used as `EATLINE(X);', this macro reads every character until it finds a
 * newline, at which points it stops; since it does not make anything of such
 * characters, it effectively advances the buffer to the next line.
 */
#define EATLINE(X) while (fgetc(X) != '\n')

#endif
