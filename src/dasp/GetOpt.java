/*
 * CDDL HEADER START
 * The contents of this file are subject to the terms
 * of the Common Development and Distribution License
 * (the License). You may not use this file except in
 * compliance with the License.
 *
 * You can obtain a copy of the License at
 * http://www.sun.com/cddl/cddl.html and legal/CDDLv1.0.txt
 * See the License for the specific language governing
 * permission and limitations under the License.
 *
 * When distributing Covered Code, include this CDDL
 * Header Notice in each file and include the License file
 * at legal/CDDLv1.0.txt.
 * If applicable, add the following below the CDDL Header,
 * with the fields enclosed by brackets [] replaced by
 * your own identifying information:
 * "Portions Copyrighted [year] [name of copyright owner]"
 *
 * Copyright 2007 Sun Microsystems Inc. All Rights Reserved
 * CDDL HEADER END
 */

package dasp;

import java.io.*;

public class GetOpt {
    
    String argList[]	= null;	    // the argument List to be parsed
    String optStr	= null;	    // the string of arguments in the form "ab:" see getopt man
    public String optArg = null;    // argument of an option
    public int optInd	= 0;	    // index of the option
    int nbArgs		= 0;	    // number of args in the argList
    int optPos		= 1;	    // position of option letter in current argument scanned
    
    static final int EOF = -1;	    // returned when no option left
    static final int ERROR = -2;    // returned on error
    
    
    public GetOpt(String[] args, String opts) {
	argList = args;
	nbArgs = argList.length;
	optStr = opts;
	optInd = 0;
	optPos = 1;
	optArg = null;
    }
    
    private void optError(String msg, char optLetter) {
	//System.err.println("Getopt: " + msg + " / " + optLetter);
    }
    
    public char getOptLetter()  {
	return argList[optInd].charAt(optPos);
    }
	
    public int getopt()  {
	optArg = null;
	
	// checking boundaries cases
	if (argList == null || optStr == null) return EOF;
	if (optInd < 0 || optInd >= nbArgs) return EOF;
	
	// get current Arg out of the argList
	String currentArg = argList[optInd];
	int argLength = currentArg.length();
	
	// checking special cases
	// if arg starts with optLetter "a", "xyz"..
	// if arg only 1 char in length "a", "-"
	if (argLength <= 1 || currentArg.charAt(0) != '-') {
	    return EOF;
	} else if (currentArg.equals("--")) { // end of options case
	    optInd++;
	    return EOF;
	}
	
	// get next "letter" from option argument
	char optLetter = currentArg.charAt(optPos);
	
	// find position of option in optStr
	int pos = optStr.indexOf(optLetter);
	if (pos == -1 || optLetter == ':') {
	    optError("illegal option", optLetter);
	    return ERROR;
	} else { // check if option requires argument
	    if (pos < optStr.length()-1 && optStr.charAt(pos+1) == ':') {
		// if remain characters after the current option in currentArg
		if (optPos != argLength-1) {
		    // take rest of current arg as the option argument
		    optArg = currentArg.substring(optPos+1);
		    optPos = argLength-1; // go to next arg below in argList
		} else { // take next arg as optArg
		    optInd++;
		    if (optInd < nbArgs &&
		    (argList[optInd].charAt(0) != '-' ||
		    argList[optInd].length() >= 2    &&
		    (optStr.indexOf(argList[optInd].charAt(1)) == -1 ||
		    argList[optInd].charAt(1) == ':')
		    )
		    ) {
			optArg = argList[optInd];
		    } else {
			optError("option '" + optLetter + "' requires an argument", optLetter);
			optArg = null;
			optLetter = '?';
		    }
		}
	    }
	}
	
	// next option argument,
	// either in currentArg or next arg in argList
	optPos++;
	
	// if no more option in currentArg
	if (optPos >= argLength) {
	    optInd++;
	    optPos = 1;  // reset postion of opt letter
	}
	return optLetter;
    }
    
}

