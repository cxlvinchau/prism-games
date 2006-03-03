//==============================================================================
//	
//	Copyright (c) 2002-2004, Andrew Hinton, Dave Parker
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

package userinterface.util;

import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import javax.swing.event.*;
import javax.swing.border.*;
import chart.*;
import userinterface.*;

/**
 *
 * @author  Andrew Hinton
 */
public class SeriesDataProperty extends SingleProperty
{
    
    private ArrayList actionListeners;
    private MultipleProperty mp;
    private boolean editingMulti;
    /** Creates a new instance of SeriesDataProperty */
    public SeriesDataProperty(PropertyOwner owner,String name, GraphList property)
    {
        this(owner, name, property, "");
    }
    
    public SeriesDataProperty(PropertyOwner owner,String name, GraphList property, String comment)
    {
        super(owner, name, property, "", false, comment);
        //FlowLayout fl = new FlowLayout(FlowLayout.CENTER, 0, 0);
        pan.setLayout(new BorderLayout());
        renderer.setBorder(null);
        pan.add(renderer);
        
        actionListeners = new ArrayList();
        
        edit = new JButton("...");
        edit.setPreferredSize(new Dimension(20, 30));
        
        //editValue = Color.BLACK;
        
        edit.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                if(editingMulti)
                {
                    if(mp != null)GraphListEditor.showEditors(GUIPrism.getGUI(), mp.getProperties());
                    fireActionPerformed(e);
                }
                else
                    getGraphList().getEditor().showEditor(null);
                fireActionPerformed(e);
            }
            
            
        });
    }
    
    
    
    public GraphList getGraphList()
    {
        return (GraphList)getProperty();
    }
    
	public void setEnabled(boolean enabled)
	{
		super.setEnabled(enabled);

		if(renderer != null)renderer.setEnabled(enabled);
		if(edit != null)edit.setEnabled(enabled);
	}
    
    
    JPanel pan = new JPanel();
    JLabel renderer = new JLabel();
    public Component getTableCellRendererComponent(JTable table, Object value,
    boolean isSelected, boolean hasFocus, int row, int column)
    {
        if(editDocked);
        {
            pan.remove(edit);
            editDocked = false;
        }
        
        renderer.setOpaque(true);
        renderer.setIcon(userinterface.GUIPrism.getIconFromImage("gridSnap.gif"));
        
        
        if (hasFocus)
        {
            pan.setBorder( UIManager.getBorder("Table.focusCellHighlightBorder") );
            if (table.isCellEditable(row, column))
            {
                pan.setForeground( UIManager.getColor("Table.focusCellForeground") );
                pan.setBackground( UIManager.getColor("Table.focusCellBackground") );
            }
        }
        else
        {
            pan.setBorder(new EmptyBorder(0, 2, 2, 1));
            pan.setForeground( UIManager.getColor("Table.focusCellForeground") );
            pan.setBackground( UIManager.getColor("Table.focusCellBackground") );
        }
        
        return pan;
    }
    
    public Component getTableCellRendererComponentMulti(JTable table, Object value,
    boolean isSelected, boolean hasFocus, int row, int column, boolean allTheSame)
    {
        if(allTheSame)
        {
            return getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
        }
        else
        {
            if(editDocked);
            {
                pan.remove(edit);
                editDocked = false;
            }
            
            renderer.setOpaque(true);
            renderer.setBackground(new Color(240,240,240));
            
            
            
            
            if (hasFocus)
            {
                pan.setBorder( UIManager.getBorder("Table.focusCellHighlightBorder") );
                if (table.isCellEditable(row, column))
                {
                    pan.setForeground( UIManager.getColor("Table.focusCellForeground") );
                    pan.setBackground( UIManager.getColor("Table.focusCellBackground") );
                }
            }
            else
            {
                pan.setBorder(new EmptyBorder(0, 2, 2, 1));
                pan.setForeground( UIManager.getColor("Table.focusCellForeground") );
                pan.setBackground( UIManager.getColor("Table.focusCellBackground") );
            }
            return pan;
        }
    }
    
    JButton edit;
    boolean editDocked = false;
    //THIS WILL NEED TO OVERRIDE THE EDITOR
    public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column)
    {
        //editValue = getColor();
        if(!editDocked)
        {
            pan.add(edit, BorderLayout.EAST);
            editDocked = true;
        }
        editingMulti = false;
        //renderer.setSelected(getBoolValue());
        pan.setBorder( UIManager.getBorder("Table.focusCellHighlightBorder") );
        
        pan.setForeground( UIManager.getColor("Table.focusCellForeground") );
        pan.setBackground( UIManager.getColor("Table.focusCellBackground") );
        
        return pan;
    }
    
    public Component getTableCellEditorComponentMulti(JTable table, Object value, boolean isSelected, int row, 
    int column, boolean allTheSame, MultipleProperty mp) 
    {
        if(allTheSame)
        {
            return getTableCellEditorComponent(table, value, isSelected, row, column);
        }
        //editValue = getColor();
        if(!editDocked)
        {
            pan.add(edit, BorderLayout.EAST);
            editDocked = true;
        }
        this.mp = mp;
        editingMulti = true;
        renderer.setBackground(new Color(240,240,240));
        //renderer.setSelected(getBoolValue());
        pan.setBorder( UIManager.getBorder("Table.focusCellHighlightBorder") );
        
        pan.setForeground( UIManager.getColor("Table.focusCellForeground") );
        pan.setBackground( UIManager.getColor("Table.focusCellBackground") );
        
        return pan;
        
    }
    
    /*public Color getEditorValue()
    {
        
        return editValue;
    }
    
    private Color editValue;
    */
    
    
    public void addListenerToEditor(ActionListener e)
    {
        actionListeners.add(e);
    }
    
    public void removeListenerFromEditor(ActionListener e)
    {
        actionListeners.remove(e);
    }
    
    public void fireActionPerformed(ActionEvent e)
    {
        for(int i = 0; i < actionListeners.size(); i++)
        {
            ((ActionListener)actionListeners.get(i)).actionPerformed(e);
        }
    }
    
    
    
}
