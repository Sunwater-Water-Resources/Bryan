
**JSON Files**
As a bit of backgournd... much of the configuration for the simulations is controlled using *JSON* files, which contain *key* and *value* pairs in a Python-like format; e.g.:
```json
{"key": "value"}
```
Where, the *value* can contain another dictionary or a list using ```[item1, item2]``` notation. The syntax in *JSON* files does not exactly mirror *Python* syntax. Main depatures include:
- The *keys* can onnly be strings - no integers or floats.
- Strings can only be bound using double inverted commas (") - no single inverted commas (') allowed.
-  Booleans use lower case for ```true``` and ```false```. 

Falty JSON syntax will often be the cause of Bryan crashes. If the error code below occurs:
```
json.decoder.JSONDecodeError: Expecting ',' delimiter:
```
This will be due to either a missing comma, a comma after the last element, or an issue with bracket closure somewhere in the JSON file. 