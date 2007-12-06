//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Andrew Hinton <ug60axh@cs.bham.ac.uk> (University of Birmingham)
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
//	
//------------------------------------------------------------------------------
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

package userinterface;

import javax.swing.*;
import java.awt.*;
import prism.*;
import java.io.*;
import userinterface.util.*;

/** This abstract class should be overridden to implement a component to be plugged
 * in to the PRISM GUI.  It extends JPanel so that it can be directly added to the
 * tabbed pane within a GUIPrism object.  A subclass should provide implementations for access
 * methods which will contain information concerning:
 * <UL>
 *	<LI>Whether the plugin should be displayed as a tab,
 *	<LI>The display text for that tab,
 *	<LI>A Menu (added to the top of the GUI)
 *	<LI>An Options Panel (added to the options dialog contained within the GUI)
 *	<LI>Any File Filters
 * </UL>
 * A subclass should also be able to handle GUIEvents by implementing the
 * processGUIEvent method, which will be called whenever a GUIEvent is generated
 * via the event handler associated with the GUIPrism object.
 * <BR>
 * This class provides access to all utilities that may be
 * required by the components:
 * <UL>
 *	<LI>Error messages
 *	<LI>Message dialogs
 *	<LI>Option dialogs
 *	<LI>File choosers
 *	<LI>Setting the taskbar text
 *	<LI>Notifying other plugins of GUIEvents
 *	<LI>Setting the state of the progress bar
 *	<LI>Making log entries
 *	<LI>Enabling this plugin
 *	<LI>Switching this plugin to the front of the gui
 *	<LI>Switching the log of the gui to the front
 * </UL>
 * @author Andrew Hinton
 * @since 21/10/03
 */
public abstract class GUIPlugin extends JPanel implements GUIEventListener, PrismSettingsListener
{
	//ATTRIBUTES
	private GUIPrism gui;
	private Prism prism;
	
	//CONSTRUCTORS
	
	/** This is the super constructor that all implementing subclasses should call, it
	 * sets up the interaction with the top level GUI, and also optionally sets up the
	 * interaction with the event
	 * handler.
	 * @param gui The parent GUI
	 * @param listens Should this plugin listen to gui events?
	 */	
	public GUIPlugin(GUIPrism gui, boolean listens)
	{
		this.gui = gui;
		this.prism = gui.getPrism();
		if(listens)gui.getEventHandler().addListener(this);
		setPreferredSize(new Dimension(800,600));
	}
	
	/** This is the super constructor that all implementing subclasses should call, it
	 * sets up the interaction with the top level GUI, and also sets up the
	 * interaction with the event
	 * handler.
	 * @param gui The parent GUI
	 */	
	public GUIPlugin(GUIPrism gui)
	{
		this(gui, true);
	}
	
	//ACCESS METHODS
	
	/** Abstract access method to be implemented to return the text that should describe
	 * this plugin in a tabbed environment.
	 * @return The tab text for this plugin.
	 */	
	public abstract String getTabText();
	
	/** Abstract access method to be implemented to provide the JMenu object that should
	 * be added to the top level GUI to provide menu functionality for this plugin.
	 * @return The JMenu to be displayed in the GUI.  Returns null if no menu is required
	 */	
	public abstract JMenu getMenu();
	
	/** Abstract access method to be implemented to provide a toolbar to be displayed in
	 * the top level GUI.
	 * @return The toolbar to be displayed in the GUI.  Returns null if no toolbar is required.
	 */	
	public abstract JToolBar getToolBar();
	
	/** Abstract access method to be implemented to provide an OptionsPanel to be
	 * displayed as part of the OptionsDialog stored in the parent GUI.
	 * @return An OptionsPanel for this plugin.  Returns null if no options panel is required.
	 */	
	public abstract OptionsPanel getOptions();
	
	/** Returns an XML representation of this plugin, for use in saving the state of the
	 * system to a file
	 * @return The XML description
	 * @deprecated This method was intended for use in projects.  Projects are no longer going to
	 * be included in the project.
	 */	
	public abstract Object getXMLSaveTree();
	
	/** Returns an XML tag for this plugin
	 * @return An XML tag
	 * @deprecated This method was intended for use in projects.  Projects are no longer going to
	 * be included in the project.
	 */	
	public abstract String getXMLIDTag();
	
	/** Returns a boolean stating whether this plugin should be displayed as a tab
	 * within the top level GUI.
	 * @return A boolean stating whether this plugin should be used as a tab.
	 */	
	public abstract boolean displaysTab();
	
	/** Utility access method to the Prism object stored within the top level user
	 * interface.
	 * @return The Prism Object.
	 */	
	public Prism getPrism()
	{
		return prism;
	}
	
	/** Returns the top level user interface
	 * @return the GUI.
	 */	
	public GUIPrism getGUI()
	{
		return gui;
	}
	
	/** Utility access method to access which component is focussed in the top level
	 * GUI.
	 * @return The focussed component in the GUI.
	 */	
	public GUIPlugin getFocussedComponent()
	{
		return gui.getFocussedPlugin();
	}
	
	/** Method to get the file selected by the filechooser
	 * @return The file from the file chooser
	 */	
	public File getChooserFile()
	{
		return gui.getChooser().getSelectedFile();
	}
	
	
	//UPDATE/UTILITY METHODS
	
	/** Do any initial actions, as rpompted by command-line args */
	public abstract void takeCLArgs(String args[]);
	
	/** When a GUIEvent is generated from some other GUIPlugin object or the parent
	 * GUIPrism object, if this plugin has been assigned to listen to events, this
	 * method is called.  It should be implemented to handle GUIEvents.  Note this
	 * method will catch GUIEvents generated from this class, so be careful to avoid
	 * recursive effects.
	 * An implementation should return false by default,
	 * unless the event should not be passed on to any more listeners,
	 * in which case it should returns true.
	 * @param e A GUIEvent from elsewhere in the userinterface.
	 */	
	public abstract boolean processGUIEvent(GUIEvent e);
	
	/** Loads an XML tree stored in 'c'.
	 * @param c The XML tree
	 * @deprecated This method was intended for use in projects.  Projects are no longer going to
	 * be included in the project.
	 */	
	public abstract void loadXML(Object c); //the object will be an xml parse tree
	
	/** A utility method which is simply a wrapper of the corresponding method in the
	 * parent GUIPrism
	 * @param title The title of the error dialog
	 * @param message The message to be displayed in the error dialog.
	 */	
	public void error(String title, String message)
	{
		gui.errorDialog(title, "Error: " + message + ".");
	}
	
	/** A utility method which is simply a wrapper of the corresponding method in the
	 * parent GUIPrism
	 * @param message The message to be displayed in the error dialog.
	 */	
	public void error(String message)
	{
		gui.errorDialog("Error: " + message + ".");
	}
	
	/** A utility method to display a message that is modal to the parent GUIPrism.
	 * @param title The string to be displayed in the title of the dialog.
	 * @param message An object containing the message to be displayed.
	 */	
	public void message(String title, Object message)
	{
		JOptionPane.showMessageDialog(gui, message, title, JOptionPane.INFORMATION_MESSAGE);
	}
	
	/** A utility method to display a message that is modal to the parent GUIPrism.
	 * @param message An object containing the message to be displayed.
	 */	
	public void message(Object message)
	{
		message("Information", message);
	}
	
	/** A utility method to display a warning message that is modal to the parent GUIPrism.
	 * @param title The string to be displayed in the title of the dialog.
	 * @param message An object containing the message to be displayed.
	 */	
	public void warning(String title, Object message)
	{
		JOptionPane.showMessageDialog(gui, "Warning: "+message+".", title, JOptionPane.WARNING_MESSAGE);
	}
	
	/** A utility method to display a warning message that is modal to the parent GUIPrism.
	 * @param message An object containing the message to be displayed.
	 */	
	public void warning(Object message)
	{
		warning("Warning", message);
	}
	
	/** A utility method to display a question dialog that is modal to the parent GUIPrism.
	 * @param title The string to be displayed in the title of the dialog.
	 * @param message An object containing the message to be displayed.
	 */	
	public int question(String title, Object message, Object[] options, int defaultOption)
	{
		if (defaultOption < 0 || defaultOption > options.length) defaultOption = 0;
		int res = JOptionPane.showOptionDialog(gui, message, title, JOptionPane.DEFAULT_OPTION , JOptionPane.QUESTION_MESSAGE, null, options, options[defaultOption]);
		if (res == JOptionPane.CLOSED_OPTION) return -1;
		return res;
	}
	
	public int question(String title, Object message, Object[] options) { return question(title, message, options, 0); } 
	public int question(Object message, Object[] options) { return question("Question", message, options, 0); }
	public int question(Object message, Object[] options, int defaultOption) { return question("Question", message, options, defaultOption); }
	
	public int questionYesNo(String title, Object message, int defaultOption)
	{
		String options[] = {"Yes", "No"};
		return question(title, message, options, defaultOption);
	}
	
	public int questionYesNo(String title, Object message) { return questionYesNo(title, message, 0); } 
	public int questionYesNo(Object message) { return questionYesNo("Question", message, 0); }
	public int questionYesNo(Object message, int defaultOption) { return questionYesNo("Question", message, defaultOption); }
	
	/** A utility method to display an option pane and return the resulting selection.
	 * @param message The message to be displayed.
	 * @param title The title of the option pane.
	 * @param buttonType The type of button to be used.  Use JOptionPane button type constants.
	 * @param messageType The type of message to be displayed. Use JOptionPane message type constants
	 * @param choices A String Array of the choices to be placed on the buttons of the option pane.
	 * @param defa The default selection string.
	 * @return Returns the index of the selected button.  The index is the corresponding index
	 * of the selection in 'choices'.
	 */	
	public int optionPane(String message, String title, int buttonType, int messageType, String[]choices, String defa)
	{
		return JOptionPane.showOptionDialog(gui, message, title, buttonType, messageType, null, choices, defa);
	}
	
	/** A utility method to show a file opening dialog with the given file filter as a
	 * default.
	 * @param ffs The list of file filters to be used within the filechooser.
	 * @param ff The file filter to be used as the default within the filechooser.
	 * @return An integer which is one of the JFileChooser selection constants.
	 */	
	public int showOpenFileDialog(GUIPrismFileFilter ffs[], GUIPrismFileFilter ff)
	{
		JFileChooser choose = gui.getChooser();
		choose.resetChoosableFileFilters();
		for(int j = 0; j < ffs.length; j++)
			choose.addChoosableFileFilter(ffs[j]);
		choose.setFileFilter(ff);
		choose.setSelectedFile(new File(""));
		return choose.showOpenDialog(gui);
	}
	
	/** A utility method to show a file saving dialog with the given file filter as a
	 * default.
	 * @param ffs The list of file filters to be used within the filechooser.
	 * @param ff The file filter to be used as the default within the filechooser.
	 * @return An integer which is one of the JFileChooser selection constants.
	 */	
	public int showSaveFileDialog(GUIPrismFileFilter ffs[], GUIPrismFileFilter ff)
	{
		JFileChooser choose = gui.getChooser();
		choose.resetChoosableFileFilters();
		for(int j = 0; j < ffs.length; j++)
			choose.addChoosableFileFilter(ffs[j]);
		choose.setFileFilter(ff);
		choose.setSelectedFile(new File(""));
		int res = choose.showSaveDialog(gui);
		if (res != JFileChooser.APPROVE_OPTION) return res;
		File file = choose.getSelectedFile();
		// check file is non-null
		if (file == null) {
			error("No file selected");
			return JFileChooser.CANCEL_OPTION;
		}
		// check for file overwrite
		if(file.exists()) {
			int selectionNo = JOptionPane.CANCEL_OPTION;
			selectionNo = optionPane("File \""+file.getPath()+"\" exists. Overwrite?", "Confirm Overwrite", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null);
			if (selectionNo != JOptionPane.OK_OPTION) return JFileChooser.CANCEL_OPTION;
		}
		return JFileChooser.APPROVE_OPTION;
	}
	
	/** A utility update method which calls the corresponding method in the parent
	 * GUIPrism object.
	 * @param message The message to be displayed on the taskbar.
	 */	
	public void setTaskBarText(String message)
	{
		gui.setTaskBarText(message);
	}
	
	/** Utility method to notify the event handler of a GUIEvent.
	 * @param e The GUIEvent to be handled.
	 */	
	public void notifyEventListeners(GUIEvent e)
	{
		gui.notifyEventListeners(e);
	}
	
	/** Utility method that calls the corresponding method in the parent GUIPrism. */	
	public void startProgress()
	{
		gui.startProgress();
	}
	
	/** Utility method that calls the corresponding method in the parent GUIPrism. */	
	public void stopProgress()
	{
		gui.stopProgress();
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 */	
	public void logln()
	{
		notifyEventListeners(new GUILogEvent(""));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(Object message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(int message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(double message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(float message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(long message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(short message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(byte message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void logln(boolean message)
	{
		notifyEventListeners(new GUILogEvent(message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(Object message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(int message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(double message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(float message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(short message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(byte message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Method to add an entry to the log contained within the parent GUIPrism.
	 * @param message The message to be added to the log
	 */	
	public void log(boolean message)
	{
		notifyEventListeners(new GUILogEvent(GUILogEvent.PRINT, message));
	}
	
	/** Utility method to change the enabled status of this GUIPlugin.
	 * @param enabled Flag to state whether this plugin should be enabled or not.
	 */	
	public void setTabEnabled(boolean enabled)
	{
		gui.enableTab(this, enabled);
	}
	
	/** Utility method to automatically jump this plugin to the front of the tabs
	 * contained within the parent GUIPrism.
	 */	
	public void tabToFront()
	{
		gui.showTab(this);
	}
	
	/** Utility method to automatically jump to the log contained within the parent
	 * GUIPrism.
	 */	
	public void logToFront()
	{
		gui.showLogTab();
	}
	
	
}
