//==============================================================================
//	
//	Copyright (c) 2005, Andrew Hinton
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package settings;
import java.awt.Color;
/**
 *
 * @author  ug60axh
 */
public abstract class ColorConstraint implements SettingConstraint
{
	
	/** Creates a new instance of ColorConstraint */
	public ColorConstraint()
	{
	}
	
	public void checkValue(Object value) throws SettingException
	{
		if(value instanceof Color)
		{
			checkValueColor((Color)value);
		}
		else
		{
			throw new SettingException("Invalid type for property, should be a Colour.");
		}
	}
	
	public abstract void checkValueColor(Color col) throws SettingException;
	
}
