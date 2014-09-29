// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "configure.hh"
#include "message.hh"
#include "string.hh"
#include "exception.hh"
#include "util.hh"
#include "types.hh"

using Util::uint;
using Util::int64;
using Util::String;
using Util::Exception;
using Util::Message;

  // So we can convert Type to strings.
  enum Type { BOOL, INT, FLOAT, STRING, ENUM };
  static const char *typeNames[] = { "bool", "int", "float", "string", "enum" };

  static const char *configKey = "config";

  union Value
  {
    bool _bool;
    int64 _int;
    double _float;
    const char *_string;
    uint _enum;
  };

  static void checkKey (const char *key);

  class Option
  {
  public:
    Option (Type type, const char *key,
            const char *description, const char *const *enumValues,
            Value value, bool printHex)
    {
      checkKey (key);
      _type = type;
      _key = strdup (key);
      _description = strdup (description);
      _enumValues = enumValues;
      _defaultValue = value;
      _value = value;
      _printHex = printHex;
    }

    bool equals (const char *key)
    {
      return (strcmp (key, _key) == 0);
    }

    Type _type;
    const char *_key;
    const char *_description;
    const char *const *_enumValues;
    Value _defaultValue;
    Value _value;
    bool _printHex;
  };

  class Map
  {
  public:
    Map ()
    {
      _count = 0;
      for (uint i = 0; i < _numBuckets; i++)
      {
        _buckets[i] = NULL;
      }
    }

    void put (Option * option)
    {
      if (get (option->_key) != NULL)
      {
        throw
          Exception (String
                     ("Configure: Option %s is already registered.",
                      option->_key));
      }
      uint index = _hash (option->_key);

      _buckets[index] = new Node (option, _buckets[index]);
      _count++;
    }

    Option *get (const char *key)
    {
      uint index = _hash (key);

      for (Node * node = _buckets[index]; node != NULL; node = node->_next)
      {
        if (node->_option->equals (key))
        {
          return node->_option;
        }
      }
      return NULL;
    }

    Option *getAlways (const char *key, Type type)
    {
      Option *option = get (key);

      if (option == NULL)
      {
        throw
          Exception (String ("Configure: Option %s is not registered.", key));
      }
      if (option->_type != type)
      {
        throw Exception (String ("Configure:: Option %s "
                                 "is not type %s, it is type %s.",
                                 key,
                                 typeNames[type], typeNames[option->_type]));
      }
      return option;
    }

    Option **getAll ()
    {
      Option **options = new Option *[_count];
      uint c = 0;

      for (uint i = 0; i < _numBuckets; i++)
      {
        for (Node * node = _buckets[i]; node != NULL; node = node->_next)
        {
          options[c++] = node->_option;
        }
      }
      assert (c == _count);
      return options;
    }

    uint count ()
    {
      return _count;
    }

  private:

    static const uint _numBuckets = 1024;

    class Node
    {

    public:
      Node (Option * option, Node * next)
      {
        _option = option;
        _next = next;
      }

      Option *_option;
      Node *_next;
    };

    uint _hash (const char *key)
    {
      uint v = 0;
      for (const char *p = key; *p != '\0'; p++)
      {
        v = ((v << 1) + v) ^ (*p);
      }
      return v % _numBuckets;
    }

    Node *_buckets[_numBuckets];
    uint _count;
  };

  // Where the option->value map lives.
  static Map map;

  // Register options...

  void
  Configure::registerBool (const char *key, bool defaultValue,
                           const char *description)
  {
    Value value;
    value._bool = defaultValue;
    Option *option = new Option (BOOL, key, description, NULL, value, false);
    map.put (option);
  }

  void
  Configure::registerInt (const char *key, int64 defaultValue,
                          const char *description, bool printHex)
  {
    Value value;
    value._int = defaultValue;
    Option *option = new Option (INT, key, description, NULL, value, printHex);
    map.put (option);
  }

  void
  Configure::registerFloat (const char *key, double defaultValue,
                            const char *description)
  {
    Value value;
    value._float = defaultValue;
    Option *option = new Option (FLOAT, key, description, NULL, value, false);
    map.put (option);
  }

  void
  Configure::registerString (const char *key, const char *defaultValue,
                             const char *description)
  {
    Value value;
    value._string = (defaultValue == NULL) ? NULL : strdup (defaultValue);
    Option *option = new Option (STRING, key, description, NULL, value, false);
    map.put (option);
  }

  void
  Configure::registerEnum (const char *key,
                           const char *const *values, uint defaultValue,
                           const char *description)
  {
    Value value;
    value._enum = defaultValue;
    Option *option = new Option (ENUM, key, description, values, value, false);
    map.put (option);
  }

  // Retrieve options...

  bool
  Configure::getBool (const char *key)
  {
    Option *option = map.getAlways (key, BOOL);
    return option->_value._bool;
  }

  int64
  Configure::getInt (const char *key)
  {
    Option *option = map.getAlways (key, INT);
    return option->_value._int;
  }

  double
  Configure::getFloat (const char *key)
  {
    Option *option = map.getAlways (key, FLOAT);
    return option->_value._float;
  }

  char *
  Configure::getString (const char *key)
  {
    Option *option = map.getAlways (key, STRING);
    if (option->_value._string == NULL)
    {
      return NULL;
    }
    else
    {
      return strdup (option->_value._string);
    }
  }

  uint
  Configure::getEnum (const char *key)
  {
    Option *option = map.getAlways (key, ENUM);

    return option->_value._enum;
  }

  // Check if the given string is valid as a key.
  static bool
  validKey (const char *key)
  {
    uint len = strlen (key);

    for (uint i = 0; i < len; i++)
    {
      int c = key[i];

      if (c == '_')
      {
        continue;
      }
      if (c == ':')
      {
        continue;
      }
      if (i == 0 && isalpha (c))
      {
        continue;
      }
      if (i > 0 && isalnum (c))
      {
        continue;
      }
      return false;
    }
    return true;
  }

  // Make sure the given string is valid as a key.
  static void
  checkKey (const char *key)
  {
    if (key == NULL)
    {
      throw Exception ("Configure: NULL key.");
    }
    if (strcmp (key, configKey) == 0)
    {
      throw
        Exception (String ("Configure: Option conflicts with builtin"
                           "config option '%s'.", configKey));
    }
    if (!validKey (key))
    {
      throw Exception (String ("Configure: Illegal key \"%s\".", key));
    }
  }

  // Return true on success, false on failure.
  // <arg> := <key>=<value>
  static bool
  parseArg (const String & arg, String & key, String & value)
  {
    // Make sure there is an '=' character in the middle of arg.
    const char *p = arg.text ();
    const char *eq = strchr (p, '=');

    if (eq == NULL)
    {
      return false;
    }
    uint offset = eq - p;

    if (offset == 0)
    {
      return false;
    }

    key.clear ();
    value.clear ();

    uint i;

    for (i = 0; i < offset; i++)
    {
      key.append (p[i]);
    }
    assert (p[i] == '=');
    for (i++; i < arg.length (); i++)
    {
      value.append (p[i]);
    }

    return true;
  }

  // Throw exception on failure.
  static void
  parseBool (const char *s, bool & val)
  {
    String str ("%s", s);

    if (str == "no")
    {
      val = false;
      return;
    }
    if (str == "off")
    {
      val = false;
      return;
    }
    if (str == "0")
    {
      val = false;
      return;
    }
    if (str == "false")
    {
      val = false;
      return;
    }

    if (str == "yes")
    {
      val = true;
      return;
    }
    if (str == "on")
    {
      val = true;
      return;
    }
    if (str == "1")
    {
      val = true;
      return;
    }
    if (str == "true")
    {
      val = true;
      return;
    }

    throw Exception (String ("Error parsing boolean value \"%s\"", s));
  }

  // Throw exception on failure.
  static void
  parseInt (const char *s, int64 & val)
  {
    int c, n;

    c = sscanf (s, "%lld%n", &val, &n);
    if (c == 1 && s[n] == '\0')
    {
      return;
    }
    c = sscanf (s, "%llx%n", &val, &n);
    if (c == 1 && s[n] == '\0')
    {
      return;
    }
    throw Exception (String ("Error parsing integer value \"%s\"", s));
  }

  // Throw exception on failure.
  static void
  parseFloat (const char *s, double &val)
  {
    int n;
    int c = sscanf (s, "%lf%n", &val, &n);

    if (c == 1 && s[n] == '\0')
    {
      return;
    }
    throw Exception (String ("Error parsing floating-point value \"%s\"", s));
  }

  // Throw exception on failure.
  static void
  parseEnum (const char *s, const char *const *values, uint & val)
  {
    String str ("%s", s);

    for (uint index = 0; values[index] != NULL; index++)
    {
      if (str == values[index])
      {
        val = index;
        return;
      }
    }
    throw Exception (String ("Error parsing enumeration value \"%s\"", s));
  }

  // Throw exception on failure.
  // Trailing whitespace is allowed.
  static Value
  parseValue (Type type, const char *const *enumValues, const char *s)
  {
    // Remove any whitespace from the end of s.
    uint len = strlen (s);

    while (len > 0 && isspace (s[len - 1]))
    {
      len--;
    }
    String val;

    val.append (len, s);

    Value value;

    switch (type)
    {
    case BOOL:
      parseBool (val, value._bool);
      break;
    case INT:
      parseInt (val, value._int);
      break;
    case FLOAT:
      parseFloat (val, value._float);
      break;
    case STRING:
      value._string = strdup (val.text ());
      break;
    case ENUM:
      parseEnum (val, enumValues, value._enum);
      break;
    default:
      assert (0);
    }
    return value;
  }

  // LEX-defined externals.
  extern int __cfl_lex ();
  extern int __cfl_lineNo;
  extern int __cfl_columnNo;
  extern FILE *__cfl_in;

  // Shared state for LEX callbaks.
  static const char *yyConfigFileName = NULL;
  static String yyKey;
  static Option *yyOption;
  enum LexStateType
  {
    lstKey,                       // Looking for a key.
    lstValue                      // Looking for a value.
  };
  static LexStateType lexState;

  // LEX callbacks.
  void
  __cfl_SyntaxError (const char *msg)
  {
    throw
      Exception (String
                 ("Syntax error while parsing configuration file:\n"
                  "%s:%d: %s", yyConfigFileName, __cfl_lineNo, msg));
  }

  void
  __cfl_FoundToken (const char *s)
  {
    switch (lexState)
    {
    case lstKey:
      {
        // Make sure the key is valid.
        if (!validKey (s))
        {
          throw
            Exception (String ("%s:%d: Invalid key \"%s\".",
                               yyConfigFileName, __cfl_lineNo, s));
        }
        // Look up the option to make sure its there.
        yyOption = map.get (s);
        if (yyOption == NULL)
        {
          Message::error(String("%s:%d: Skipping unknown option \"%s\".\n",
                                 yyConfigFileName, __cfl_lineNo, s));
  //              throw Exception (String (
  //                  "%s:%d: Unknown option \"%s\".",
  //                  yyConfigFileName, __cfl_lineNo, s));
        }
        yyKey.clear ();
        yyKey.append ("%s", s);
        lexState = lstValue;
        break;
      }
    case lstValue:
      {
        // Parse the value and add it to the map.
        if (yyOption != NULL)
        {
          try
          {
            yyOption->_value =
              parseValue (yyOption->_type, yyOption->_enumValues, s);
          }
          catch (Exception & e)
          {
            throw
              Exception (String ("%s:%d: %s for \"%s\".",
                                 yyConfigFileName, __cfl_lineNo,
                                 e.msg (), yyKey.text ()));
          }
        }
        lexState = lstKey;
        break;
      }
    default:
      assert (0);
    }
  }
  void
  __cfl_FoundString (const char *s)
  {
    __cfl_FoundToken (s);
  }

  void
  __cfl_AllDone ()
  {
    switch (lexState)
    {
    case lstKey:
      {
        // Ok.
        break;
      }
    case lstValue:
      {
        throw
          Exception (String ("%s:%d: Missing value for \"%s\".",
                             yyConfigFileName, __cfl_lineNo, yyKey.text ()));
        break;
      }
    default:
      assert (0);
    }
  }

  // Parse the configuration file using a LEX-generated lexer.
  static void
  processConfigFile (const char *fileName)
  {
    FILE *strm = NULL;

    try
    {
      strm = Util::openInputStrm (fileName);
      yyConfigFileName = fileName;
      lexState = lstKey;
      __cfl_in = strm;
      int rval = __cfl_lex ();

      if (rval != 0)
      {
        throw Exception (String ("Lex error (%d) on file %s.", rval, fileName));
      }
    }
    catch (Exception & e)
    {
      if (strm != NULL)
      {
        fclose (strm);
      }
      throw;
    }
  }

  void
  Configure::init (const char *key, const char *value)
  {
    assert (key != NULL);

    // Look for special config key.
    if (strcmp (key, configKey) == 0)
    {
      processConfigFile (value);
      return;
    }
    // Look up the key.
    Option *option = map.get (key);

    if (option == NULL)
    {
      throw Exception (String ("Option %s not found.", key));
    }
    // Parse the value and store the result.
    try
    {
      option->_value = parseValue (option->_type, option->_enumValues, value);
    }
    catch (Exception & e)
    {
      throw Exception (String ("%s for option \"%s\".", e.msg (), key));
    }
  }

  // Apply all the options on the command-line to the map.
  // Read any config files specified on the command-line as well.
  uint
  Configure::init (int argc, const char **argv)
  {
    // First scan command-line for the 'Config' flags.
    for (int index = 1; index < argc; index++)
    {
      String arg ("%s", argv[index]);
      String key, value;

      if (arg == "--")
      {
        index++;
        break;
      }
      if (parseArg (arg, key, value))
      {
        if (!validKey (key))
        {
          break;
        }
      }
      else if (argv[index][0] == '-')
      {
        key.append ("%s", argv[index] + 1);
        if (!validKey (key))
        {
          break;
        }
        if (index == argc - 1)
        {
          break;
        }
        index++;
        value.append ("%s", argv[index]);
      }
      else
      {
        break;
      }
      if (key == configKey)
      {
        processConfigFile (value);
      }
    }

    // Process options on the command-line.
    int index;

    for (index = 1; index < argc; index++)
    {
      String arg ("%s", argv[index]);
      String key, value;
      Option *option;

      if (arg == "--")
      {
        index++;
        break;
      }
      if (parseArg (arg, key, value))
      {
        assert (validKey (key));
        if (key == configKey)
        {
          continue;
        }
        option = map.get (key);
        if (option == NULL)
        {
          break;
        }
      }
      else if (argv[index][0] == '-')
      {
        key.append ("%s", argv[index] + 1);
        if (!validKey (key))
        {
          break;
        }
        assert (validKey (key));
        if (index == argc - 1)
        {
          break;
        }
        if (key == configKey)
        {
          index++;
          continue;
        }
        option = map.get (key);
        if (option == NULL)
        {
          break;
        }
        index++;
        value.append ("%s", argv[index]);
      }
      else
      {
        break;
      }
      try
      {
        option->_value = parseValue (option->_type, option->_enumValues,
                                     value.text ());
      }
      catch (Exception & e)
      {
        throw
          Exception (String ("%s for option \"%s\".", e.msg (), key.text ()));
      }
    }

    // Return the index at which command-line processing failed.
    // Return argc on no error.
    return index;
  }

  // For qsort().
  static int
  cmpOption (const void *x, const void *y)
  {
    Option *a = *(Option **) x;
    Option *b = *(Option **) y;
    bool aRoot = (strstr (a->_key, "::") == NULL);
    bool bRoot = (strstr (a->_key, "::") == NULL);

    if (aRoot && !bRoot)
    {
      return 1;
    }
    if (!aRoot && bRoot)
    {
      return -1;
    }
    return strcmp (a->_key, b->_key);
  }

  // Append the given value to the string.
  static void
  valueString (String & val, Option * option, Value value)
  {
    switch (option->_type)
    {
    case BOOL:
      val.append ("%s", value._bool ? "true" : "false");
      break;
    case INT:
      if (option->_printHex)
      {
        val.append ("0x%llx", value._int);
      }
      else
      {
        val.append ("%lld", value._int);
      }
      break;
    case FLOAT:
      val.append ("%g", value._float);
      break;
    case STRING:
      val.append ("%s", (value._string == NULL) ? "<null>" : value._string);
      break;
    case ENUM:
      val.append ("%s", option->_enumValues[value._enum]);
      break;
    default:
      assert (0);
    }
  }

  // A little state machine to do word wrapping.
#define CHAROUT(c) { \
      assert ((c) != '\n'); \
      fputc ((c), strm); \
      col++; \
  }
#define NEWLINE { \
      fputc ('\n', strm); \
      col = 1; \
      for (uint i = 0; i < indent - 1; i++) { \
          CHAROUT(' '); \
      } \
  }
  static void
  wrapText (FILE * strm, uint col, uint indent, uint width, const char *text)
  {
    const char *start = text;

    while (*start != '\0')
    {
      // Start a new line if we hit a newline character.
      if (*start == '\n')
      {
        NEWLINE;
        start++;
        continue;
      }
      // For whitespace, echo spaces as long as they fit on the current line.
      if (isspace (*start))
      {
        if (col < width)
        {
          CHAROUT (' ');
        }
        start++;
        continue;
      }
      // Find the end of the next word.
      const char *end = start;

      while (*end != '\0' && !isspace (*end))
      {
        end++;
      }
      uint len = end - start;

      // Output the word if it fits on this line.
      if (col + len < width)
      {
        while (start != end)
        {
          CHAROUT (*start);
          start++;
        }
        continue;
      }
      // The word goes off the end of the line.
      // If we're in the middle of the line, then put this word on the next line.
      if (col > indent)
      {
        NEWLINE;
        continue;
      }
      // The word is wider than the output column.  We need to break the word
      // in half.
      while (col < width)
      {
        CHAROUT (*start);
        start++;
        assert (start <= end);
      }
    }
    fputc ('\n', strm);
  }

  // Print a usage message to strm.
  void
  Configure::usage()
  {
    FILE* strm = stderr;
    uint argWidth = 25;
    String fmt ("  %%-%ds ", argWidth);
    uint indent = 4 + argWidth;
    uint width = 79;

    fprintf (strm, "\n");
    String key ("%s[=<null>]", configKey);

    fprintf (strm, fmt.text (), key.text ());
    wrapText (strm, indent, indent, width, "Configuration file.");

    Option **options = map.getAll ();

    qsort ((void *) options, map.count (), sizeof (*options), cmpOption);
    for (uint i = 0; i < map.count (); i++)
    {
      assert (options[i] != NULL);
      String val;

      valueString (val, options[i], options[i]->_defaultValue);
      String key ("%s[=%s]", options[i]->_key, val.text ());
      String desc ("%s", options[i]->_description);

      if (options[i]->_type == ENUM)
      {
        String values ("%s", "Possible values are {");

        for (uint index = 0; options[i]->_enumValues[index] != NULL; index++)
        {
          if (index > 0)
          {
            values.append (',');
          }
          values.append ("%s", options[i]->_enumValues[index]);
        }
        values.append ("}.");
        desc.append ("  %s", values.text ());
      }
      uint col = (key.length () > argWidth) ? width : indent;

      fprintf (strm, fmt.text (), key.text ());
      wrapText (strm, col, indent, width, desc.text ());
    }
    delete[]options;
  }

  // Print the map to strm.
  void
  Configure::show()
  {
    Message::debug("\n\nCONFIGURATION BEGIN",1);
    Option **options = map.getAll ();
    qsort ((void *) options, map.count (), sizeof (*options), cmpOption);
    for (uint i = 0; i < map.count (); i++)
    {
      assert (options[i] != NULL);
      String val;

      valueString (val, options[i], options[i]->_value);
      String key ("%s", options[i]->_key);

      Message::debug(String("%s %s", options[i]->_key, val.text ()),1);
    }
    Message::debug("CONFIGURATION END\n\n",1);
    delete[]options;
  }


