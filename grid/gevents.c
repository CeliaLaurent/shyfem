
/************************************************************************\
 *
 *    Copyright (C) 1994-1995,1998  Georg Umgiesser
 *
 *    This file is part of SHYFEM.
 *
 *    SHYFEM is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SHYFEM is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with SHYFEM. Please see the file COPYING in the main directory.
 *    If not, see <http://www.gnu.org/licenses/>.
 *
 *    Contributions to this file can be found below in the revision log.
 *
\************************************************************************/


/************************************************************************\
 *
 * gevents.c - event routines for gcc under DOS
 *
 * revision log :
 *
 * 04.12.1994	ggu	routines written from scratch
 * 10.12.1994	ggu	routines adapted to gcc through Event Queue Library
 * 06.12.1995	ggu	TranslateKeyboardEvent for special keys
 * 13.02.1998	ggu	adjourned to GRX2.0
 *
\************************************************************************/


#include <stdio.h>

#define GRX_VERSION_GGU	2

#if GRX_VERSION_GGU == 1
#include <grx.h>
#include <mousex.h>
#else
#include <grx20.h>
#endif

#include <keys.h>

#include "general.h"
#include "mouse.h"
#include "keybd.h"

#include "events.h"
/*
#include "eventque.h"
*/

#if GRX_VERSION_GGU == 1
#define	GrMouseInit	MouseInit
#define	GrMouseUnInit	MouseUnInit
#define	GrMouseGetEvent(E,e) 	MouseGetEvent( (E) , (e) )
#endif

#define GR_M_BUTTON_LEFT		( GR_M_LEFT_DOWN | GR_M_LEFT_UP )
#define GR_M_BUTTON_RIGHT		( GR_M_RIGHT_DOWN | GR_M_RIGHT_UP )
#define GR_M_BUTTON_MIDDLE		( GR_M_MIDDLE_DOWN | GR_M_MIDDLE_UP )

static QEventMask DefaultEventMask	= QKeyPressMask
					| QButtonPressMask
					| QExposureMask
					;
static QEventMask ActualEventMask;

/*
static EventQueue *ActualEventQueue;
static EventRecord ActualEvent;
*/
static GrMouseEvent ActualEvent;

static int TranslateKeyboardEvent( int key );

#define GGU_DEBUG 1
#undef GGU_DEBUG
#ifdef GGU_DEBUG
#include "debug.h"
#endif

/*********************************************************************/

void QInitEvent( void )

{
	QSelectEvent( DefaultEventMask );
/*	ActualEventQueue = EventQueueInit(200,0,0); */
	GrMouseInit();
/*	MouseShow(); */
}

void QDeInitEvent( void )

{
	GrMouseUnInit();
}

QEventMask QSelectedEvent( void )

{
	return ActualEventMask;
}

void QSelectEvent( QEventMask eventmask )

{
	ActualEventMask = eventmask;
}

void QAddEvent( QEventMask eventmask )

{
	ActualEventMask = ActualEventMask | eventmask;
}

void QDeleteEvent( QEventMask eventmask )

{
	ActualEventMask = ActualEventMask & ~eventmask;
}

void QNextEvent( QEvent *eventp )

{
	unsigned char button;
	unsigned char moved;
	int event;
	static int first = TRUE;

	if( first ) {

		first = FALSE;
		eventp->type = QExpose;

		return;
	}

	eventp->type = QNullEvent;

/*
	moved = ActualEventQueue->evq_moved;
	event = EventQueueNextEvent( ActualEventQueue , &ActualEvent );
*/
	GrMouseGetEvent( GR_M_EVENT , &ActualEvent );

	moved = ActualEvent.flags & GR_M_MOTION;
	event = ActualEvent.flags & ( GR_M_KEYPRESS | GR_M_BUTTON_CHANGE );
	button = ActualEvent.flags & GR_M_BUTTON_CHANGE;

	if( !event && moved && ( ActualEventMask & QPointerMoveMask ) ) {

		/* only if no other event is present */

		eventp->type = QPointerMove;
		eventp->button.x = ActualEvent.x;
		eventp->button.y = ActualEvent.y;

		return;
	}

	if( !event ) return;

	if( 	      button
		 && ( ActualEventMask & QButtonPressMask ) 
		 ) {

		eventp->type = QButtonPress;
		eventp->button.x = ActualEvent.x;
		eventp->button.y = ActualEvent.y;


		if( button & GR_M_BUTTON_LEFT ) {

			eventp->button.button = QButtonLeft;
			if( button & GR_M_BUTTON_DOWN )
				eventp->button.press = QButtonDown; 
			else
				eventp->button.press = QButtonUp; 

		} else if( button & GR_M_BUTTON_RIGHT ) {

			eventp->button.button = QButtonRight;
			if( button & GR_M_BUTTON_DOWN )
				eventp->button.press = QButtonDown; 
			else
				eventp->button.press = QButtonUp; 

		} else {

			eventp->button.button = QButtonMiddle;
			if( button & GR_M_BUTTON_DOWN )
				eventp->button.press = QButtonDown; 
			else
				eventp->button.press = QButtonUp; 

		}

#ifdef GGU_DEBUG
	GDS("Button Event : "); 
	GDI(button);
	GDI(eventp->button.x);
	GDI(eventp->button.y);
	GDNL();
#endif

		return;
	}

	if( 	    ( event & GR_M_KEYPRESS )
		 && ( ActualEventMask & QKeyPressMask ) 
		 ) {

		eventp->type = QKeyPress;
		eventp->key  = TranslateKeyboardEvent(ActualEvent.key);

#ifdef GGU_DEBUG
	GDS("Keyboard Event : "); GDI(eventp->key); GDNL();
#endif

		return;
	}
}
		
static int TranslateKeyboardEvent( int key )

{
	if( key >= 0x20 && key <= 0x7e ) return key;

	if( key == K_BackSpace ) 
		return QKeyBackSpace;
	if( key == K_Tab ) 
		return QKeyTab;
	if( key == K_Return ) 
		return QKeyReturn;
	if( key == K_Escape ) 
		return QKeyEscape;

	if( key == K_EUp ) 
		return QKeyUp;
	if( key == K_EDown ) 
		return QKeyDown;
	if( key == K_ELeft ) 
		return QKeyLeft;
	if( key == K_ERight ) 
		return QKeyRight;

	return QKeyUnknown;
}
