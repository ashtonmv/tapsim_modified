/************************************************************************
*                                                                       *
* TAPSim - an atom probe data simulation program                        *
*                                                                       *
* Copyright (C) 2011 Christian Oberdorfer                               *
*                                                                       *
* This program is free software: you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation, either version 3 of the License, or any *
* any later version.                                                    *
*                                                                       *
* This program is distributed in the hope that it will be useful,       *
* but WITHOUT ANY WARRANTY; without even the implied warranty of        *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
* GNU General Public License for more details.                          *
*                                                                       *
* You should have received a copy of the GNU General Public License     *
* along with this program.  If not, see 'http://www.gnu.org/licenses'   *
*                                                                       *
************************************************************************/

#include "file_util.h"

std::string File_Utils::trim(const std::string& source)
{
	const char* whiteValues = " \t\n";

	const size_t start = source.find_first_not_of(whiteValues,0);
	if (start == std::string::npos) return std::string();
	
	const size_t end = source.find_last_not_of(whiteValues,std::string::npos);
	if (end == std::string::npos) return std::string();
	
	return source.substr(start,end-start+1);
}
